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
        self.output_files = []  # Store the output filenames for multiple models

        for model in self.cov_model_path:
            output_file = f"output_cmscan_{os.path.basename(model)}.txt"
            cmd = f"cmscan -o {output_file} {model} input_for_cm_scan.fa"
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            process.communicate()
            self.output_files.append(output_file)  # Store the output file name

    def _parse_result_file(self):
        def _parse_chunk(chunk):
            list_cmscan_hits = []
            flag_in = False
            for line in chunk:
                if "----   --------- ------" in line:
                    flag_in = True
                    continue
                if "inclusion threshold" in line or "No hits detected that satisfy reporting thresholds":
                    flag_in = False

                if "Hit alignments:" in line:
                    break

                if flag_in:
                    line_elements = line.split()
                    if line_elements:
                        start = line_elements[6]
                        end = line_elements[7]
                        strand = line_elements[8]
                        e_val = float(line_elements[2])
                        score = float(line_elements[3])
                        complete_hit = CMScanHit(start, end, strand, e_val, score)
                        list_cmscan_hits.append(complete_hit)

            return list_cmscan_hits

        self.list_hits_best = []  # Reset hits before processing multiple files

        for output_file in self.output_files:  # Loop through all output files
            if not os.path.exists(output_file):
                continue  # Skip missing files

            with open(output_file) as f:
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

                self.list_hits_best.append(best_hit)  # Store the best hit per sequence

    def _clean_up(self):
        try:
            os.remove("input_for_cm_scan.fa")  # Remove the input file

            # Remove all cmscan output files
            for output_file in self.output_files:
                if os.path.exists(output_file):
                    os.remove(output_file)

        except Exception as e:
            print(f"Cleanup error: {e}")  # Print error if cleanup fails

    def output(self):
        return self.list_hits_best