import os
import subprocess
from os.path import join, basename
import shutil
from modules.consistency_score_maker import get_header_column_index



from modules.consistency_score_maker import get_header_column_index

import os
from os.path import join, basename

from modules.consistency_score_maker import get_header_column_index


class CMScanRunCompleteDNA:
    def __init__(self, dna_path, model_path, folder_intermediate_files):
        self.dna_path = dna_path
        self.model_path = model_path
        self.folder_intermediate_files = folder_intermediate_files

        self.dna = None
        self.acc_num = basename(dna_path).split(".")[0]

        self.intervals_scores = []

        self.merged_intervals_forward = []
        self.merged_intervals_reversed = []

        self._make_intermediate_files_folder()
        self._get_dna()
        self._run_cm_scan()
        self._parse_result_file()
        self._clean_up()

    @staticmethod
    def reverse_com(seq):
        dict_complementary = {"A": "T", "G": "C", "T": "A", "C": "G", "N": "N"}
        return "".join([dict_complementary.get(x, x) for x in seq])[::-1]

    @staticmethod
    def format_e_value(e_value):
        return f"{e_value:.2e}"  # Convert to scientific notation

    def _make_intermediate_files_folder(self):
        os.makedirs(self.folder_intermediate_files, exist_ok=True)

    def _get_dna(self):
        with open(self.dna_path) as f:
            self.dna = "".join(line.strip() for line in f if not line.startswith(">"))

    def _run_cm_scan(self):
        for model in self.model_path:
            output_file = f"output_cmscan_{basename(model)}.txt"
            cmd = ["cmscan", "-o", output_file, model, self.dna_path]
            subprocess.run(cmd, check=True)

    def _parse_result_file(self):
        self.intervals_scores = []
        for model in self.model_path:
            output_file = f"output_cmscan_{basename(model)}.txt"
            if not os.path.exists(output_file):
                continue
            with open(output_file) as f:
                lines = f.readlines()

            flag_in = False
            for line in lines:
                if "----   --------- ------" in line:
                    flag_in = True
                    continue
                if "inclusion threshold" in line or "No hits detected" in line:
                    flag_in = False
                if "Hit alignments:" in line:
                    break
                if flag_in:
                    line_elements = line.split()
                    if line_elements:
                        start, end = int(line_elements[6]), int(line_elements[7])
                        strand, e_val = line_elements[8], float(line_elements[2])
                        if start > end:
                            start, end = end, start
                            strand = "-" if strand == "+" else "+"
                        self.intervals_scores.append((start, end, strand, e_val))

        self.intervals_scores.sort()

    def _clean_up(self):
        for model in self.model_path:
            output_file = f"output_cmscan_{basename(model)}.txt"
            if os.path.exists(output_file):
                os.remove(output_file)

    def write_raw_hits_to_csv(self, csv_file_name):
        with open(csv_file_name, "w") as f:
            f.write("acc_num,start,end,e_value,hit_sequence\n")
            for start, end, strand, e_value in self.intervals_scores:
                seq = self.dna[start - 1:end] if strand == "+" else self.reverse_com(self.dna[start - 1:end])
                f.write(f"{self.acc_num},{start},{end},{self.format_e_value(e_value)},{seq}\n")

    def compute_and_write_merged_predictions(self, csv_file_name):
        def merge_intervals(intervals):
            intervals.sort()
            merged = []
            for interval in intervals:
                start, end, strand, e_val = interval
                if merged and start - merged[-1][1] <= 150 and strand == merged[-1][2]:
                    prev_start, prev_end, prev_strand, prev_eval = merged.pop()
                    merged.append((min(prev_start, start), max(prev_end, end), prev_strand, min(prev_eval, e_val)))
                else:
                    merged.append(interval)
            return merged

        forward_intervals = [(s, e, strand, e_val) for s, e, strand, e_val in self.intervals_scores if strand == "+"]
        reverse_intervals = [(s, e, strand, e_val) for s, e, strand, e_val in self.intervals_scores if strand == "-"]

        self.merged_intervals_forward = merge_intervals(forward_intervals)
        self.merged_intervals_reversed = merge_intervals(reverse_intervals)

        with open(csv_file_name, "w") as f:
            f.write("acc_num,start,end,best_e_value,hit_sequence\n")
            for start, end, strand, best_e_value in self.merged_intervals_forward:
                seq = self.dna[start - 1:end]
                f.write(f"{self.acc_num},{start},{end},{self.format_e_value(best_e_value)},{seq}\n")
            for start, end, strand, best_e_value in self.merged_intervals_reversed:
                seq = self.reverse_com(self.dna[start - 1:end])
                f.write(f"{self.acc_num},{start},{end},{self.format_e_value(best_e_value)},{seq}\n")



def filter_csv_file_model_run(csv_file_name, output_file_name, e_value_threshold, hit_length_threshold):
    with open(csv_file_name, "r") as f:
        lines = f.readlines()
        header = lines[0]
        index_acc_num = get_header_column_index(header, ",", "acc_num")
        index_start = get_header_column_index(header, ",", "start")
        index_end = get_header_column_index(header, ",", "end")
        index_hit_sequence = get_header_column_index(header, ",", "hit_sequence")
        index_e_value = get_header_column_index(header, ",", "e_value")
        lines = lines[1:]
        dict_lines = {}
        filtered_lines = []
        for line in lines:
            line = line.strip()
            line_elements = line.split(",")
            acc_num = line_elements[index_acc_num]
            start = line_elements[index_start]
            end = line_elements[index_end]
            hit_sequence = line_elements[index_hit_sequence]
            e_value = float(line_elements[index_e_value])

            if e_value > e_value_threshold:
                filtered_lines.append(line)
                continue
            if len(hit_sequence) < hit_length_threshold:
                filtered_lines.append(line)
                continue

            if acc_num not in dict_lines:
                dict_lines[acc_num] = {(start, end): (e_value, line)}
            else:
                if (start, end) not in dict_lines[acc_num]:
                    dict_lines[acc_num][(start, end)] = (e_value, line)
                else:
                    if e_value < dict_lines[acc_num][(start, end)][0]:
                        filtered_lines.append(dict_lines[acc_num][(start, end)][1])
                        dict_lines[acc_num][(start, end)] = (e_value, line)
                    else:
                        filtered_lines.append(line)

    with open(output_file_name, "w") as f:
        f.write(header)
        for acc_num in dict_lines:
            for interval in dict_lines[acc_num]:
                f.write(dict_lines[acc_num][interval][1])
                f.write("\n")
        f.write("\n\n Filtered out candidates (duplicates; short hits; high e-value): \n")
        for line in filtered_lines:
            f.write(line)
            f.write("\n")

