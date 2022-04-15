import os


class IntarnaRun:
    def __init__(self, repeat_sequence, tracr_sequence):
        self.repeat_sequence = repeat_sequence
        self.tracr_sequence = tracr_sequence

        self._run_intarna()
        self._parse_intarna_result()
        self._clean_up()

    def _run_intarna(self):
        with open("repeat.fa", "w") as f:
            f.write(f">repeat\n{self.repeat_sequence}\n")

        with open("tracr.fa", "w") as f:
            f.write(f">tracr\n{self.tracr_sequence}\n")

        cmd = "IntaRNA -t tracr.fa -q repeat.fa --out intarna_result.csv --outMode C"
        os.system(cmd)

    def _parse_intarna_result(self):
        with open("intarna_result.csv") as f:
            lines = f.readlines()
            if len(lines) > 1:
                info_line = lines[1]

                line_elements = info_line.strip().split(";")
                repeat_interval = "_".join([line_elements[4], line_elements[5]])
                self.tracr_interval = "_".join([line_elements[1], line_elements[2]])
                self.interaction = "____".join([line_elements[6], line_elements[7]])
                self.score = float(line_elements[-1])
            else:
                self.tracr_interval = ""
                self.interaction = ""
                self.score = ""

    def _clean_up(self):
        try:
            os.remove("repeat.fa")
            os.remove("tracr.fa")
            os.remove("intarna_result.csv")
        except Exception:
            pass

    def output(self):
        return self.tracr_interval, self.interaction, self.score