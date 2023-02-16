from os import listdir
from os.path import isfile, join

from modules.crispr_identify_summary_parser import CRISPRIdentifySummaryParser
from modules.anti_repeat_search import AntiRepeatSearch, AntiRepeatFilter
from modules.candidates_based_on_anti_repeat import TracrRNASecondPartSearch
from modules.intarna_run import IntarnaRun
from modules.complete_terminator_search import TerminatorCompleteSearch
from modules.hmm_model_run import CMScanOnCandidates
from modules.general_helpers import reverse_com
from modules.consistency_score_maker import consistency_score_maker
from modules.run_identify_and_identifyer import run_crispr_identify, run_crispr_cas_identifier_on_folder
from modules.general_helpers import get_closest_interval, string_tuple_to_tuple
from modules.candidate_ranking import candidate_ranking
from modules.folder_maker import folder_maker
from modules.model_run_only import CMScanRunCompleteDNA
from modules.pattern_search import anti_repeat_search_type_v
from modules.model_run_only import filter_csv_file_model_run


class CompleteTracrSearchWithModel:
    def __init__(self, folder_output, folder_genome,
                 temp_folder_path, output_file_name,
                 blast_similarity, blast_coverage, hmm_model, dict_weights):

        self.folder_output = folder_output
        self.folder_genome = folder_genome
        self.temp_folder_path = temp_folder_path
        self.output_file_name = output_file_name
        self.blast_similarity = blast_similarity
        self.blast_coverage = blast_coverage
        self.hmm_model = hmm_model
        self.dict_weights = dict_weights

        self._pipeline_run()

    def _pipeline_run(self):
        print("Creating folders if necessary...")
        folder_maker(self.folder_output, self.temp_folder_path)

        print("Running CRISPR array search, it may take a while...")
        folder_result_identify = "ResultsIdentify"
        run_crispr_identify(self.folder_genome, folder_result_identify)
        crispr_identify_summary_file = join(folder_result_identify, "Complete_summary.csv")
        dict_arrays = CRISPRIdentifySummaryParser(crispr_identify_summary_file).output()

        print("Running CRISPR cas search, it may take a while...")
        dict_all_cas = run_crispr_cas_identifier_on_folder(self.folder_genome)
        dict_all_cas = {key.split(".")[0]: value for key, value in dict_all_cas.items()}

        files_folder = [f for f in listdir(self.folder_genome) if isfile(join(self.folder_genome, f))]
        print(f"\t\tWorking with {len(files_folder)} files \n\n")

        for index, file_name in enumerate(files_folder, 1):
            print(f"\t\tWorking with file {index} - {file_name}")
            full_path_genome = join(self.folder_genome, file_name)

            acc_num = file_name.split(".")[0]
            if acc_num not in dict_arrays:
                continue
            arrays = dict_arrays[acc_num]
            consensus_repeats = [cr.consensus for cr in arrays]
            consensus_repeats_alternative_orientation = [reverse_com(r) for r in consensus_repeats]
            array_intervals = [(cr.start, cr.end) for cr in arrays]

            print(f"1. Searching for anti repeats")
            ar_search = AntiRepeatSearch(consensus_repeats, full_path_genome, self.temp_folder_path)
            list_fasta_hits = ar_search.output()

            ar_search_alternative_orientation = AntiRepeatSearch(consensus_repeats_alternative_orientation,
                                                                 full_path_genome, self.temp_folder_path)

            list_fasta_hits_alternative_orientation = ar_search_alternative_orientation.output()

            print("2. Filtering Anti-repeats")
            ar_filter = AntiRepeatFilter(path_dna_file=full_path_genome,
                                         list_array_intervals=array_intervals,
                                         list_consensus_repeats=consensus_repeats,
                                         list_hits=list_fasta_hits,
                                         repeat_coverage=self.blast_coverage,
                                         blast_similarity=self.blast_similarity,
                                         folder_intermediate_files=self.temp_folder_path)

            ar_filter_alternative_orientation = AntiRepeatFilter(path_dna_file=full_path_genome,
                                                                 list_array_intervals=array_intervals,
                                                                 list_consensus_repeats=consensus_repeats_alternative_orientation,
                                                                 list_hits=list_fasta_hits_alternative_orientation,
                                                                 repeat_coverage=self.blast_coverage,
                                                                 blast_similarity=self.blast_similarity,
                                                                 folder_intermediate_files=self.temp_folder_path)

            anti_repeat_candidates, anti_repeat_candidates_wrong = ar_filter.output()
            anti_repeat_candidates_alt, anti_repeat_candidates_alt_wrong = ar_filter_alternative_orientation.output()

            print("3. Forming initial TracrRNA candidates")
            anti_all_ways = [anti_repeat_candidates, anti_repeat_candidates_wrong,
                             anti_repeat_candidates_alt, anti_repeat_candidates_alt_wrong]

            tracr_candidates_all_ways = []

            for anti in anti_all_ways:
                tracr_rna_search = TracrRNASecondPartSearch(anti, full_path_genome)
                tracr_candidates = tracr_rna_search.output()
                tracr_candidates_all_ways.append(tracr_candidates)

            print("4. Improving the results using CRISPRRNA-TracrRNA interaction prediction")
            intarna_result_all_ways = []
            for tracr_candidates in tracr_candidates_all_ways:
                list_intarna_results = []
                for tracr_candidate in tracr_candidates:
                    repeat_sequence = tracr_candidate.anti_repeat.crispr_array_consensus_repeat
                    tracr_sequence = tracr_candidate.complete_seq
                    ir = IntarnaRun(repeat_sequence, tracr_sequence)
                    tracr_interval, interaction, score = ir.output()
                    list_intarna_results.append((tracr_interval, interaction, score))

                intarna_result_all_ways.append(list_intarna_results)

            print("5. Searching for the terminator sequences")
            terminators_all_ways = []
            for tracr_candidates in tracr_candidates_all_ways:
                all_tracr_sequences = [tr_candidate.complete_seq for tr_candidate in tracr_candidates]
                terminator_cs = TerminatorCompleteSearch(all_tracr_sequences)
                output_dictionary_by_column = terminator_cs.output_by_column()
                output_all_terminators_by_column = terminator_cs.output_by_columns_all()
                terminator_location = output_dictionary_by_column["terminator location"]
                terminator_best_score = output_dictionary_by_column["terminator score"]
                terminator_all_locations = output_all_terminators_by_column["all terminator locations"]
                terminator_all_scores = output_all_terminators_by_column["all terminator scores"]

                terminator_all_info = [terminator_all_locations, terminator_all_scores,
                                       terminator_location, terminator_best_score]

                terminators_all_ways.append(terminator_all_info)

            print("6. Searching the hits with the provided tail model")
            all_scan_results = []
            for tracr_candidates in tracr_candidates_all_ways:
                all_tracr_sequences_global_window = [tr_candidate.anti_repeat_upstream + tr_candidate.anti_repeat.anti_repeat_seq + tracr_candidate.downstream_region
                                                     for tr_candidate in tracr_candidates]

                cm_scan_on_candidates = CMScanOnCandidates(all_tracr_sequences_global_window, self.hmm_model)
                scan_results = cm_scan_on_candidates.output()
                all_scan_results += scan_results

            print("7. Writing results")
            path_to_output = join(self.folder_output, acc_num + ".csv")
            with open(path_to_output, "w") as f:
                header = ",".join(["accession_number",
                                   "crispr_array_index",
                                   "crispr_array_category",
                                   "crispr_array_score",
                                   "crispr_array_start",
                                   "crispr_array_end",
                                   "crispr_array_repeat_consensus",
                                   "crispr_array_orientation",
                                   "crispr_orientation_flag",
                                   "anti_repeat_sequence",
                                   "anti_repeat_start",
                                   "anti_repeat_end",
                                   "anti_repeat_direction",
                                   "anti_repeat_relative_location",
                                   "anti_repeat_distance_from_crispr_array",
                                   "anti_repeat_similarity",
                                   "anti_repeat_coverage",
                                   "anti_repeat_similarity_coverage_multiplication",
                                   "anti_repeat_upstream",
                                   "tracr_rna_taken_flag",
                                   "tracr_rna_tail_sequence",
                                   "tracr_rna_global_window_sequence",
                                   "tracr_rna_sequence",
                                   "intarna_anti_repeat_interaction_interval",
                                   "intarna_anti_repeat_interaction",
                                   "poli_u_signal_coordinates",
                                   "interaction_energy",
                                   "terminator_all_locations",
                                   "terminator_all_scores",
                                   "best_terminator_location",
                                   "best_terminator_score",
                                   "terminator_presence_flag",
                                   "tail_model_hit_location",
                                   "tail_model_hit_score",
                                   "tail_presence_flag",
                                   "closest_corresponding_cas_interval",
                                   "distance_to_cas"])

                f.write(f"{header}\n")

                for way_index in range(4):
                    array_orientation = "Predicted" if (way_index in [0, 1]) else "Alternative"
                    tracr_rna_flag_taken = "5'-3'" if (way_index in [0, 2]) else "3'-5'"
                    tracr_candidates = tracr_candidates_all_ways[way_index]
                    intarna_results = intarna_result_all_ways[way_index]
                    terminator_results = terminators_all_ways[way_index]

                    pair = [tracr_candidates, intarna_results]
                    for tracr_index, (tracr_candidate, intarna_result) in enumerate(zip(*pair)):
                        anti_repeat_candidate = tracr_candidate.anti_repeat
                        crispr_array_index = anti_repeat_candidate.crispr_array_origin_index
                        crispr_array_info = dict_arrays[acc_num][crispr_array_index]
                        accession_number = acc_num
                        crispr_array_index = crispr_array_index
                        crispr_array_category = crispr_array_info.category
                        crispr_array_score = crispr_array_info.score
                        crispr_array_start = crispr_array_info.start
                        crispr_array_end = crispr_array_info.end
                        crispr_array_repeat_consensus = anti_repeat_candidate.crispr_array_consensus_repeat
                        crispr_array_orientation_flag = array_orientation
                        if array_orientation == "Predicted":
                            crispr_array_orientation = crispr_array_info.strand
                        else:
                            crispr_array_original_orientation = crispr_array_info.strand
                            crispr_array_orientation = list({"Forward", "Reversed"} -
                                                            {crispr_array_original_orientation})[0]
                        anti_repeat_sequence = anti_repeat_candidate.anti_repeat_seq
                        anti_repeat_start = anti_repeat_candidate.anti_repeat_start
                        anti_repeat_end = anti_repeat_candidate.anti_repeat_end
                        anti_repeat_direction = anti_repeat_candidate.anti_repeat_strand
                        anti_repeat_relative_location = anti_repeat_candidate.anti_repeat_relative_location
                        anti_repeat_similarity = anti_repeat_candidate.fasta_similarity
                        anti_repeat_coverage = anti_repeat_candidate.fasta_coverage
                        anti_repeat_similarity_coverage_multiplication = anti_repeat_similarity * anti_repeat_coverage
                        anti_repeat_distance_from_crispr_array = anti_repeat_candidate.anti_repeat_distance_from_crispr_array
                        anti_repeat_upstream = tracr_candidate.anti_repeat_upstream
                        tracr_rna_tail_sequence = tracr_candidate.downstream_region
                        tracr_rna_global_window_sequence = anti_repeat_upstream + anti_repeat_sequence + tracr_rna_tail_sequence
                        tracr_rna_sequence = anti_repeat_sequence + tracr_rna_tail_sequence
                        intarna_anti_repeat_interaction_interval = intarna_result[0]
                        intarna_anti_repeat_interaction = intarna_result[1]
                        poli_u_tail = str(tracr_candidate.signal_coordinates).replace(",", "-")
                        interaction_energy = intarna_result[2]
                        terminator_all_locations = str(terminator_results[0][tracr_index]).replace(",", "-")
                        terminator_all_scores = str(terminator_results[1][tracr_index]).replace(",", "-")
                        best_terminator_location = str(terminator_results[2][tracr_index]).replace(",", "-")
                        best_terminator_score = str(terminator_results[3][tracr_index]).replace(",", "-")

                        terminator_presence_flag = 1 if terminator_results[0][tracr_index] != "NA" else 0

                        scan_result = all_scan_results[tracr_index]
                        if scan_result:
                            scan_interval = f"{scan_result.Start}_{scan_result.End}"
                            scan_score = str(scan_result.Score)
                        else:
                            scan_interval = "NA"
                            scan_score = "NA"

                        tail_presence_flag = 1 if scan_result else 0

                        if acc_num not in dict_all_cas:
                            all_cas_intervals = []
                        else:
                            all_cas_intervals = [key for key, value in dict_all_cas[acc_num].items() if value == "cas9" or value == "cas12"]
                        #all_cas_intervals = [string_tuple_to_tuple(interval) for interval in all_cas_intervals]

                        interval_tracr_rna_anti_repeat = (int(anti_repeat_start), int(anti_repeat_end))
                        cas_distance, closest_cas_interval = get_closest_interval(interval_tracr_rna_anti_repeat, all_cas_intervals)
                        closest_cas_interval = str(closest_cas_interval).replace(", ", "-")

                        complete_line = ",".join([str(x) for x in [accession_number,
                                                                   crispr_array_index,
                                                                   crispr_array_category,
                                                                   crispr_array_score,
                                                                   crispr_array_start,
                                                                   crispr_array_end,
                                                                   crispr_array_repeat_consensus,
                                                                   crispr_array_orientation,
                                                                   crispr_array_orientation_flag,
                                                                   anti_repeat_sequence,
                                                                   anti_repeat_start,
                                                                   anti_repeat_end,
                                                                   anti_repeat_direction,
                                                                   anti_repeat_relative_location,
                                                                   anti_repeat_distance_from_crispr_array,
                                                                   anti_repeat_similarity,
                                                                   anti_repeat_coverage,
                                                                   anti_repeat_similarity_coverage_multiplication,
                                                                   anti_repeat_upstream,
                                                                   tracr_rna_flag_taken,
                                                                   tracr_rna_tail_sequence,
                                                                   tracr_rna_global_window_sequence,
                                                                   tracr_rna_sequence,
                                                                   intarna_anti_repeat_interaction_interval,
                                                                   intarna_anti_repeat_interaction,
                                                                   poli_u_tail,
                                                                   interaction_energy,
                                                                   terminator_all_locations,
                                                                   terminator_all_scores,
                                                                   best_terminator_location,
                                                                   best_terminator_score,
                                                                   terminator_presence_flag,
                                                                   scan_interval,
                                                                   scan_score,
                                                                   tail_presence_flag,
                                                                   closest_cas_interval,
                                                                   cas_distance]])

                        f.write(f"{complete_line}\n")

        print("\n\t\tCreating the final_summary")
        header = ""
        final_content = []
        all_summary_files = [f for f in listdir(self.folder_output) if isfile(join(self.folder_output, f))]
        for index, s_file in enumerate(all_summary_files):
            with open(join(self.folder_output, s_file)) as fr:
                lines = fr.readlines()
                if index == 0:
                    header = lines[0]
                final_content += lines[1:]

        with open(self.output_file_name, "w") as f:
            f.write(header)
            for line in final_content:
                f.write(line)

        print("\n\t\tAdding consistency scores to the final summary file")
        consistency_score_maker(self.output_file_name, ",", self.output_file_name)

        print("\n\t\tRanking the final summary file")
        candidate_ranking(self.output_file_name, ",", self.dict_weights, self.output_file_name)


class TracrSearchWihtModelOnly:
    def __init__(self, folder_input, folder_output, summary_file_name, path_to_model, flag_perform_anti_repeat_search):
        self.folder_input = folder_input
        self.folder_output = folder_output
        self.summary_file_name = summary_file_name
        self.path_to_model = path_to_model
        self.flag_perform_anti_repeat_search = flag_perform_anti_repeat_search

        self._run_model_search()

    def _run_model_search(self):
        folder_maker(self.folder_output, "")
        print("\n\t\tRunning the model search")
        all_files = [f for f in listdir(self.folder_input) if isfile(join(self.folder_input, f))]
        for index, file_name in enumerate(all_files):
            print(f"\tWorking with file {file_name}  {index+1} out of {len(all_files)}")
            full_path_to_file = join(self.folder_input, file_name)
            cm_scan_c_dna = CMScanRunCompleteDNA(full_path_to_file, self.path_to_model,  self.folder_output)
            cm_scan_c_dna.report_csv_output_file(join(self.folder_output, file_name + "_report.csv"))


        lines_report = []
        files_in_temp_folder = [f for f in listdir(self.folder_output)
                                if isfile(join(self.folder_output, f))]

        for index, file_name in enumerate(files_in_temp_folder):
            with open(join(self.folder_output, file_name)) as f:
                if index == 0:
                    lines_report = f.readlines()
                else:
                    lines_report += f.readlines()[1:]

        with open(self.summary_file_name, "w") as fw:
            fw.writelines(lines_report)

        if self.flag_perform_anti_repeat_search:
            print("\n\t\tPerforming anti-repeat search")
            anti_repeat_search_type_v(self.summary_file_name, ",", self.summary_file_name)

        filter_csv_file_model_run(self.summary_file_name, self.summary_file_name, 0.0001, 30)
