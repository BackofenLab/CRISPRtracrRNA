import os
import subprocess

from collections import namedtuple
CMScanHit = namedtuple('CMScanHit', ["Start", "End", "Strand", "Eval", "Score"])


class CMScanOnCandidates:
    def __init__(self, list_sequences, cov_model_path):
        self.list_sequences = list_sequences
        self.cov_model_path = cov_model_path

        self.list_hits_all = []
        self.list_hits_best = []

        self._create_input_file()
        self._run_cm_scan()
        self._parse_result_file()
        self._clean_up()

    def _create_input_file(self):
        with open("input_for_cm_scan.fa", "w") as f:
            for seq_index, sequence in enumerate(self.list_sequences):
                header = f">{seq_index}\n"
                f.write(header)
                f.write(sequence)
                f.write("\n")

    def _run_cm_scan(self):
        cmd = f"cmscan -o output_cmscan.txt {self.cov_model_path} input_for_cm_scan.fa"
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        a, b = process.communicate()

    def _parse_result_file(self):
        def _parse_chunk(chunk):
            list_cmscan_hits = []
            flag_in = False
            for line in chunk:
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
                        e_val = float(line_elements[2])
                        score = float(line_elements[3])
                        complete_hit = CMScanHit(start, end, strand, e_val, score)
                        list_cmscan_hits.append(complete_hit)

            return list_cmscan_hits

        with open("output_cmscan.txt") as f:
            lines = f.readlines()

        indexes_headers = [index for index, line in enumerate(lines) if "Query:" in line]
        chunks = [lines[index_s:index_e] for index_s, index_e in zip(indexes_headers, indexes_headers[1:])]
        last_chunk = lines[indexes_headers[-1]:]
        chunks.append(last_chunk)

        for index, chunk in enumerate(chunks):
            list_model_hits = _parse_chunk(chunk)
            if list_model_hits:
                best_hit = max(list_model_hits, key=lambda x: x.Score)
            else:
                best_hit = None

            self.list_hits_best.append(best_hit)

    def _clean_up(self):
        try:
            os.remove("input_for_cm_scan.fa")
            os.remove("output_cmscan.txt")
        except Exception:
            pass

    def output(self):
        return self.list_hits_best


if __name__ == "__main__":
    list_seqs = ["AAGGCUUUGUCCGUACACAACUUGAAAAAGUGCGCACCGAUUCGGUGCUUUUUU", "AAGGCUUUGUCCGUACACAACUUGAAAAAGUGCGCACCGAUUCGGUGCUUUUU", "AAGGCAGUGAUUUUUAAUCCAGUCCGUAUUCAGCUUGAAAAAGUGAGCACCGAUUCGGUGCUUUUUU", "ACCCAGAAAATA"]
    model = "/home/alex/PycharmProjects/CRISPRtracr/tools/cov_matrix/s2_b14/cov_matrix.cm"
    results = CMScanOnCandidates(list_seqs, model).output()
    print(results)