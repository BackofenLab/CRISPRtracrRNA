import os
from os.path import isfile, join
import shutil
from collections import namedtuple


AntiRepeat = namedtuple('AntiRepeat', ["crispr_array_origin_index",
                                       "crispr_array_consensus_repeat",
                                       "anti_repeat_seq",
                                       "anti_repeat_start",
                                       "anti_repeat_end",
                                       "anti_repeat_strand",
                                       "anti_repeat_relative_location",
                                       "anti_repeat_distance_from_crispr_array",
                                       "correct_wrong_orientation",
                                       "fasta_similarity",
                                       "fasta_coverage"])


class AntiRepeatSearch:
    def __init__(self, list_consensus_repeats,
                 path_dna_file, folder_intermediate_files,):
        self.list_consensus_repeats = list_consensus_repeats
        self.path_dna_file = path_dna_file
        self.folder_intermediate_files = folder_intermediate_files

        self.fasta_hits = []

        self._create_input_file()
        self._run_fasta()
        self._parse_fasta()
        self._clean_up()

    def _create_input_file(self):
        with open("all_repeat_candidates.fa", "w") as f:
            for index, consensus_repeat in enumerate(self.list_consensus_repeats):
                f.write(f">{index}\n")
                f.write(f"{consensus_repeat}\n")

    def _run_fasta(self):
        cmd = f"tools/fasta-36.3.8g/bin/fasta36 {self.path_dna_file} " \
              f"all_repeat_candidates.fa  -m 8 > fasta_similarity.fastab"
        os.system(cmd)
        cmd = f"tools/fasta-36.3.8g/bin/fasta36 all_repeat_candidates.fa {self.path_dna_file} " \
              f" -m 8 > fasta_similarity_inverted.fastab"
        os.system(cmd)

    def _parse_fasta(self):
        with open("fasta_similarity.fastab") as f:
            self.list_fasta_hits = [line.strip().split() for line in f.readlines()]

        with open("fasta_similarity_inverted.fastab") as f:
            list_fasta_hits_inverted = [line.strip().split() for line in f.readlines()]

        list_fasta_hits_inverted_correct_val_order = []
        for h in list_fasta_hits_inverted:
            if int(h[6]) < int(h[7]):
                hit_new = [h[1], h[0], h[2], h[3], h[4], h[5], h[8], h[9], h[6], h[7], h[10], h[11]]
            else:
                hit_new = [h[1], h[0], h[2], h[3], h[4], h[5], h[9], h[8], h[7], h[6], h[10], h[11]]
            list_fasta_hits_inverted_correct_val_order.append(hit_new)

        self.list_fasta_hits += list_fasta_hits_inverted_correct_val_order

    def _clean_up(self):
        try:
            shutil.move('all_repeat_candidates.fa', join(self.folder_intermediate_files, "all_repeat_candidates.fa"))
            shutil.move('fasta_similarity.fastab', join(self.folder_intermediate_files, "fasta_similarity.fastab"))
            shutil.move('fasta_similarity_inverted.fastab', join(self.folder_intermediate_files,
                                                                 "fasta_similarity_inverted.fastab"))
        except Exception:
            pass

    def output(self):
        return self.list_fasta_hits


class AntiRepeatFilter:
    def __init__(self, path_dna_file, list_array_intervals, list_consensus_repeats,
                 list_hits, repeat_coverage, blast_similarity,
                 folder_intermediate_files):

        self.path_dna_file = path_dna_file
        self.list_array_intervals = list_array_intervals
        self.list_consensus_repeats = list_consensus_repeats
        self.list_fasta_hits = list_hits
        self.repeat_coverage = repeat_coverage
        self.blast_similarity = blast_similarity
        self.folder_intermediate_files = folder_intermediate_files

        self.anti_repeats_correct = []
        self.anti_repeats_wrong = []

        self._read_dna_file()
        self._filter_out_hits_by_coverage_and_similarity()
        self._filter_out_hits_by_arrays()
        self._filter_out_duplicates()
        self._fill_in_dictionary_of_containers()

    def _read_dna_file(self):
        with(open(self.path_dna_file)) as f:
            lines = f.readlines()
        self.dna = "".join([line.strip() for line in lines[1:]])

    def _filter_out_hits_by_coverage_and_similarity(self):
        hits_filtered_by_coverage_and_similarity = []

        for fasta_hit in self.list_fasta_hits:
            hit_length = int(fasta_hit[3])
            blast_sim_score = float(fasta_hit[2]) * 0.01
            header = int(fasta_hit[1])
            corresponding_repeat = self.list_consensus_repeats[header]
            repeat_length = len(corresponding_repeat)

            if hit_length > (repeat_length * self.repeat_coverage):
                if blast_sim_score > self.blast_similarity:
                    hits_filtered_by_coverage_and_similarity.append(fasta_hit)

        self.list_fasta_hits = hits_filtered_by_coverage_and_similarity

    def _filter_out_hits_by_arrays(self):
        list_fasta_hits_filtered_by_arrays = []
        for fasta_hit in self.list_fasta_hits:
            hit_start = int(fasta_hit[6])
            for interval in self.list_array_intervals:
                start, end = interval
                if start <= hit_start <= end:
                    break
            else:
                list_fasta_hits_filtered_by_arrays.append(fasta_hit)

        self.list_fasta_hits = list_fasta_hits_filtered_by_arrays

    def _filter_out_duplicates(self):
        hit_signatures = []
        filtered_duplications_hits = []
        for hit in self.list_fasta_hits:
            hit_signature = f"{hit[1]}_{hit[6]}_{hit[7]}"
            if hit_signature not in hit_signatures:
                hit_signatures.append(hit_signature)
                filtered_duplications_hits.append(hit)

        self.list_fasta_hits = filtered_duplications_hits

    def _write_intermediate_anti_repeat_hits(self):
        with open(join(self.folder_intermediate_files, "fasta_hits_after_filtration.txt"), "w") as f:
            for hit in self.list_fasta_hits:
                line = "\t".join(hit)
                f.write(line)
                f.write("\n")

    @staticmethod
    def reverse_com(seq):
        dict_complementary = {"A": "T", "G": "C", "T": "A", "C": "G", "N": "N"}
        new_seq = "".join([dict_complementary.get(x, x) for x in seq])
        new_seq = new_seq[::-1]
        return new_seq

    def _fill_in_dictionary_of_containers(self):
        for fasta_hit in self.list_fasta_hits:
            crispr_array_origin_index = int(fasta_hit[1])
            crispr_array_consensus_repeat = self.list_consensus_repeats[crispr_array_origin_index]
            array_start, array_end = self.list_array_intervals[crispr_array_origin_index]

            hit_fasta_start = int(fasta_hit[6])
            hit_fasta_end = int(fasta_hit[7])

            q_hit_start = int(fasta_hit[8])
            q_hit_end = int(fasta_hit[9])

            fasta_similarity = float(fasta_hit[2]) * 0.01
            fasta_coverage = int(fasta_hit[3]) / len(crispr_array_consensus_repeat)

            hit_start, hit_end = sorted([hit_fasta_start, hit_fasta_end])

            if (hit_fasta_start < hit_fasta_end) and (q_hit_start < q_hit_end):
                flag_hit_direction = "reversed"
                flag_different_strand_hit_direction = "forward"
            else:
                flag_hit_direction = "forward"
                flag_different_strand_hit_direction = "reversed"

            if flag_hit_direction == "forward":
                hit_sequence = self.dna[hit_start - 1:hit_end - 1]

                hit_different_strand = self.reverse_com(hit_sequence)
            else:
                hit_sequence = self.dna[hit_start - 1:hit_end - 1]
                hit_sequence = self.reverse_com(hit_sequence)

                hit_different_strand = self.dna[hit_start - 1:hit_end - 1]

            flag_relative_direction = "Upstream" if hit_start < array_start else "Downstream"
            if flag_relative_direction == "Upstream":
                anti_repeat_distance_from_crispr_array = array_start - hit_end
            else:
                anti_repeat_distance_from_crispr_array = hit_start - array_end

            an = AntiRepeat(crispr_array_origin_index=crispr_array_origin_index,
                            crispr_array_consensus_repeat=crispr_array_consensus_repeat,
                            anti_repeat_seq=hit_sequence,
                            anti_repeat_start=hit_start,
                            anti_repeat_end=hit_end,
                            anti_repeat_strand=flag_hit_direction,
                            anti_repeat_relative_location=flag_relative_direction,
                            anti_repeat_distance_from_crispr_array=anti_repeat_distance_from_crispr_array,
                            correct_wrong_orientation="correct",
                            fasta_similarity=fasta_similarity,
                            fasta_coverage=fasta_coverage)

            an_wrong = AntiRepeat(crispr_array_origin_index=crispr_array_origin_index,
                                  crispr_array_consensus_repeat=crispr_array_consensus_repeat,
                                  anti_repeat_seq=hit_different_strand,
                                  anti_repeat_start=hit_start,
                                  anti_repeat_end=hit_end,
                                  anti_repeat_strand=flag_different_strand_hit_direction,
                                  anti_repeat_relative_location=flag_relative_direction,
                                  anti_repeat_distance_from_crispr_array=anti_repeat_distance_from_crispr_array,
                                  correct_wrong_orientation="wrong",
                                  fasta_similarity=fasta_similarity,
                                  fasta_coverage=fasta_coverage)

            self.anti_repeats_correct.append(an)
            self.anti_repeats_wrong.append(an_wrong)

    def output(self):
        return self.anti_repeats_correct, self.anti_repeats_wrong
