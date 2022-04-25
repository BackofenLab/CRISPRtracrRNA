import subprocess
import os
from collections import namedtuple

CRISPRArray = namedtuple('CRISPRArray', ["start",
                                         "end",
                                         "consensus",
                                         "strand",
                                         "score",
                                         "category"])


class CRISPRIdentifySummaryParser:
    def __init__(self, summary_file_path):
        self.summary_file = summary_file_path
        self.dict_array_info = {}

        self._parse_file()

    def _parse_file(self):
        with open(self.summary_file) as f:
            lines = f.readlines()
        for line in lines[1:]:
            line_elements = line.strip().split(",")
            acc_number = line_elements[0]
            start = int(line_elements[4])
            end = int(line_elements[5])
            consensus = line_elements[7]
            strand = line_elements[11]
            score = float(line_elements[13])
            category = line_elements[12]
            if category in ["Bona-fide", "Possible"]:
                cr_array = CRISPRArray(start, end, consensus, strand, score, category)
                if acc_number in self.dict_array_info:
                    self.dict_array_info[acc_number].append(cr_array)
                else:
                    self.dict_array_info[acc_number] = [cr_array]

    def output(self):
        return self.dict_array_info
