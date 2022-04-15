import subprocess
import os


class TerminatorCompleteSearch:
    def __init__(self, list_of_tracr_rna_sequences):
        self.list_of_tracr_rna_sequences = list_of_tracr_rna_sequences

        self._write_into_file()
        self._run_erpin()
        self._parse_result()
        self._clean_up()
        self._reformat_as_columns()

    def _write_into_file(self):
        with open("temp_erpin_file.fa", "w") as f:
            for index, value in enumerate(self.list_of_tracr_rna_sequences):
                f.write(f">{index}\n")
                f.write(f"{value}\n")

    def _run_erpin(self):
        exe_path = "tools/erpin/erpin"
        database_path = "tools/erpin/rho-indep.epn"
        input_file_path = "temp_erpin_file.fa"
        output_result = "tools/erpin/erpin_result.txt"
        cmd = f"{exe_path} {database_path}  {input_file_path} -1,4 -add 2 4 2 -pcw 3.0 -cutoff 100%% > {output_result}"

        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        a, b = process.communicate()

    def _parse_result(self):
        with open("tools/erpin/erpin_result.txt") as f:
            lines = f.readlines()

        indexes_headers = [index for index, line in enumerate(lines) if line[0] == ">"]
        headers = [int(lines[index].strip()[1:]) for index in indexes_headers]
        dict_best_hits = {}
        dict_all_hits = {}
        for header_index, header in zip(indexes_headers, headers):
            index_hit_line = header_index + 1
            hit_line_description = lines[index_hit_line]
            hit_line_elements = hit_line_description.strip().split()
            hit_line_strand = hit_line_elements[0]
            if hit_line_strand == "FW":
                start_end = hit_line_elements[2]
                hit_start = int(start_end.split("..")[0])
                hit_end = int(start_end.split("..")[1])
                score = float(hit_line_elements[4])
                if header not in dict_best_hits:
                    dict_best_hits[header] = (hit_start, hit_end), score
                    dict_all_hits[header] = [((hit_start, hit_end), score)]
                else:
                    dict_all_hits[header].append(((hit_start, hit_end), score))
                    existing_hit_end, existing_score = dict_best_hits[header]
                    if score < existing_score:
                        dict_best_hits[header] = (hit_start, hit_end), score

        self.dict_all_hits = dict_all_hits
        self.dict_best_hits = dict_best_hits

    def _reformat_as_columns(self):
        self.dict_column_output = {}
        self.dict_by_column_all = {}
        list_location = []
        list_score = []

        for line_index in range(len(self.list_of_tracr_rna_sequences)):
            if line_index not in self.dict_best_hits:
                list_location.append("NA")
                list_score.append("NA")
            else:
                list_location.append(self.dict_best_hits[line_index][0])
                list_score.append(self.dict_best_hits[line_index][1])

        self.dict_column_output["terminator location"] = list_location
        self.dict_column_output["terminator score"] = list_score

        list_all_locations = []
        list_all_scores = []

        for line_index in range(len(self.list_of_tracr_rna_sequences)):
            if line_index not in self.dict_best_hits:
                list_all_locations.append("NA")
                list_all_scores.append("NA")
            else:
                list_all_locations.append([x[0] for x in self.dict_all_hits[line_index]])
                list_all_scores.append([x[1] for x in self.dict_all_hits[line_index]])

        self.dict_by_column_all["all terminator locations"] = list_all_locations
        self.dict_by_column_all["all terminator scores"] = list_all_scores

    def _clean_up(self):
        try:
            os.remove("temp_erpin_file.fa")
        except Exception:
            pass

        try:
            os.remove("tools/erpin/erpin_result.txt")
        except Exception:
            pass

    def output(self):
        return self.dict_all_hits, self.dict_best_hits

    def output_by_column(self):
        return self.dict_column_output

    def output_by_columns_all(self):
        return self.dict_by_column_all

