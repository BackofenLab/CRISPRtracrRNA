import os
from os.path import join
import shutil
from os.path import basename
from modules.consistency_score_maker import get_header_column_index


class CMScanRunCompleteDNA:
    def __init__(self, dna_path, model_path, folder_intermediate_files):
        self.dna_path = dna_path
        self.model_path = model_path
        self.folder_intermediate_files = folder_intermediate_files

        self.dna = None
        self.acc_num = basename(dna_path).split(".")[0]

        self.intervals_scores = []

        self.scores_forward = []
        self.scores_reversed = []

        self.merged_intervals_forward = []
        self.merged_intervals_reversed = []

        self.merged_intervals_sequences_forward = []
        self.merged_intervals_sequences_reversed = []

        self._make_intermediate_files_folder()
        self._get_dna()
        self._run_cm_scan()
        self._parse_result_file()
        self._clean_up()
        self._compute_overlaps()
        self._get_sequences_for_intervals()

    @staticmethod
    def reverse_com(seq):
        dict_complementary = {"A": "T", "G": "C", "T": "A", "C": "G", "N": "N"}
        new_seq = "".join([dict_complementary.get(x, x) for x in seq])
        new_seq = new_seq[::-1]
        return new_seq

    def _make_intermediate_files_folder(self):
        try:
            os.mkdir(self.folder_intermediate_files)
        except OSError:
            pass

    def _get_dna(self):
        with open(self.dna_path) as f:
            lines = f.readlines()

        self.dna = "".join(line.strip() for line in lines[1:])

    def _run_cm_scan(self):
        cmd = f"tools/infernal/binaries/cmscan -o output_cmscan.txt {self.model_path} {self.dna_path}"
        os.system(cmd)

    def _parse_result_file(self):
        with open("output_cmscan.txt") as f:
            lines = f.readlines()

        flag_in = False
        for line in lines:
            if "----   --------- ------" in line:
                flag_in = True
                continue
            if "inclusion threshold" in line:
                flag_in = False

            if "No hits detected that satisfy reporting thresholds" in line:
                flag_in = False

            if "Hit alignments:" in line:
                break

            if flag_in:
                line_elements = line.split()
                if line_elements:
                    line_elements = line.split()
                    start = line_elements[6]
                    end = line_elements[7]
                    strand = line_elements[8]
                    e_val = line_elements[2]
                    score = line_elements[3]
                    complete_hit = f"{start}_{end}_{strand}_{e_val}_{score}"
                    self.intervals_scores.append(complete_hit)

    def _clean_up(self):
        try:
            #shutil.move('output_cmscan.txt', join(self.folder_intermediate_files, "output_cmscan.txt"))
            os.remove("output_cmscan.txt")
        except Exception:
            pass

    def _compute_overlaps(self):
        def merging_two_intervals(interval_one, interval_two):
            start_one, end_one = interval_one
            start_two, end_two = interval_two

            if (start_two <= end_one <= end_two) or (start_two <= start_one <= end_two):
                start = min(start_one, start_two)
                end = max(end_one, end_two)
                return start, end

        def merging_intervals(intervals, evals):
            for interval_one in intervals:
                for interval_two in intervals:
                    if interval_one != interval_two:
                        merged_interval = merging_two_intervals(interval_one, interval_two)
                        if merged_interval:
                            index_interval1 = intervals.index(interval_one)
                            intervals.remove(interval_one)
                            e_val1 = evals[index_interval1]
                            evals.pop(index_interval1)

                            index_interval2 = intervals.index(interval_two)
                            intervals.remove(interval_two)
                            e_val2 = evals[index_interval2]
                            evals.pop(index_interval2)

                            intervals.append(merged_interval)
                            evals.append(min(e_val1, e_val2))

                            return False, intervals, evals

            return True, intervals, evals

        def iterative_merging_intervals(intervals, e_vals):
            flag_merged_all_possible = False

            while not flag_merged_all_possible:
                flag_merged_all_possible, intervals, e_vals = merging_intervals(intervals, e_vals)

            return intervals, e_vals

        cmscan_e_vals_forward = [float(str_interval.split("_")[3]) for str_interval in self.intervals_scores if
                                 str_interval.split("_")[2] == "+"]
        cmscan_intervals_forward = [[int(str_interval.split("_")[0]), int(str_interval.split("_")[1])]
                                    for str_interval in self.intervals_scores if str_interval.split("_")[2] == "+"]

        cmscan_intervals_forward = [[min(interval), max(interval)] for interval in cmscan_intervals_forward]

        cmscan_intervals_forward_merged, cmscan_evals_forward_merged = iterative_merging_intervals(cmscan_intervals_forward,
                                                                                                   cmscan_e_vals_forward)
        self.merged_intervals_forward = cmscan_intervals_forward_merged

        cmscan_e_vals_reversed = [float(str_interval.split("_")[3]) for str_interval in self.intervals_scores if
                                  str_interval.split("_")[2] == "-"]
        cmscan_intervals_reversed = [[int(str_interval.split("_")[0]), int(str_interval.split("_")[1])]
                                     for str_interval in self.intervals_scores if str_interval.split("_")[2] == "-"]

        cmscan_intervals_reversed = [[min(interval), max(interval)] for interval in cmscan_intervals_reversed]
        cmscan_intervals_reversed_merged, cmscan_evals_reversed_merged = iterative_merging_intervals(cmscan_intervals_reversed,
                                                                                                     cmscan_e_vals_reversed)

        self.merged_intervals_forward = cmscan_intervals_forward_merged
        self.merged_intervals_reversed = cmscan_intervals_reversed_merged

        self.merged_intervals_e_vals_forward = cmscan_evals_forward_merged
        self.merged_intervals_e_vals_reversed = cmscan_evals_reversed_merged

    def _get_sequences_for_intervals(self):
        for interval in self.merged_intervals_forward:
            start, end = interval
            extended_start = max(0, start-30)
            extended_end = min(end+30, len(self.dna))
            seq = self.dna[start-1:end-1]
            extended_seq = self.dna[extended_start:extended_end]
            self.merged_intervals_sequences_forward.append((interval, seq, extended_seq))

        for interval in self.merged_intervals_reversed:
            start, end = interval
            extended_start = max(0, start - 30)
            extended_end = min(end + 30, len(self.dna))
            seq = self.dna[start - 1:end - 1]
            seq = self.reverse_com(seq)
            extended_seq = self.dna[extended_start:extended_end]
            extended_seq = self.reverse_com(extended_seq)
            self.merged_intervals_sequences_reversed.append((interval, seq, extended_seq))

    def output_intervals(self):
        return self.intervals_scores

    def output_merged_intervals(self):
        return self.merged_intervals_forward, self.merged_intervals_reversed

    def output_merged_intervals_sequences(self):
        return self.merged_intervals_sequences_forward, self.merged_intervals_sequences_reversed

    def report_cm_scan_results(self):
        file_path = join(self.folder_intermediate_files, "cmscan_report.txt")
        with open(file_path, "w") as f:
            for interval_seq in self.merged_intervals_sequences_forward:
                interval, seq, extended_seq = interval_seq
                start, end = interval
                header = f">{start}_{end}_forward\n"
                f.write(header)
                f.write(seq)
                f.write("\n")

            for interval_seq in self.merged_intervals_sequences_reversed:
                interval, seq, extended_seq = interval_seq
                start, end = interval
                header = f">{start}_{end}_reversed\n"
                f.write(header)
                f.write(seq)
                f.write("\n")

    def report_csv_output_file(self, csv_file_name):
        with open(csv_file_name, "w") as f:
            header = ",".join(["acc_num", "start", "end", "hit_sequence", "e_value", "extended_sequence"]) + "\n"
            f.write(header)

            for interval_seq, e_value in zip(self.merged_intervals_sequences_forward, self.merged_intervals_e_vals_forward):
                interval, seq, extended_seq = interval_seq
                start, end = interval
                line = ",".join([str(x) for x in [self.acc_num, start, end, seq, e_value, extended_seq]]) + '\n'
                f.write(line)

            for interval_seq, e_value in zip(self.merged_intervals_sequences_reversed, self.merged_intervals_e_vals_reversed):
                interval, seq, extended_seq = interval_seq
                start, end = interval
                line = ",".join([str(x) for x in [self.acc_num, start, end, seq, e_value, extended_seq]]) + '\n'
                f.write(line)


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

