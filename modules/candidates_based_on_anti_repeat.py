import re
from collections import namedtuple

AntiRepeatTracrRNA = namedtuple("AntiRepeatTracrRNA", ["anti_repeat",
                                                       "complete_seq",
                                                       "anti_repeat_upstream",
                                                       "anti_repeat_seq",
                                                       "downstream_region",
                                                       "anti_repeat_upstream_coordinates",
                                                       "anti_repeat_start",
                                                       "anti_repeat_end",
                                                       "downstream_region_coordinates",
                                                       "signal_coordinates",
                                                       "anti_repeat_direction"])


class TracrRNASecondPartSearch:
    def __init__(self, list_anti_repeats, path_dna_file):
        self.list_anti_repeats = list_anti_repeats
        self.path_dna_file = path_dna_file

        self.list_tracr_rna = []

        self._read_fasta_file()
        self._search_all_tracr_rna()

    def _read_fasta_file(self):
        with(open(self.path_dna_file)) as f:
            lines = f.readlines()
        self.dna = "".join([line.strip() for line in lines[1:]])

    @staticmethod
    def reverse_com(seq):
        dict_complementary = {"A": "T", "G": "C", "T": "A", "C": "G", "N": "N"}
        new_seq = "".join([dict_complementary.get(x, x) for x in seq])
        new_seq = new_seq[::-1]
        return new_seq

    @staticmethod
    def get_subseqs_by_signal(sequence, characters):
        all_tuples_match_coordinates = []
        for character in characters:
            pattern = character + "{4,}"
            matched_strings = re.findall(pattern, sequence)
            re_iter = re.finditer(pattern, sequence)
            start_end = [(m.start(0), m.end(0)) for m in re_iter]
            tuples_match_coordinates = list(zip(matched_strings, start_end))
            all_tuples_match_coordinates += tuples_match_coordinates

        sorted_tuples = sorted(all_tuples_match_coordinates, key=lambda x: (len(x[0]) - (x[1][0]) / 10000),
                               reverse=True)
        return sorted_tuples

    def _search_all_tracr_rna(self):
        for anti_repeat in self.list_anti_repeats:
            anti_repeat_seq = anti_repeat.anti_repeat_seq
            hit_start = anti_repeat.anti_repeat_start
            hit_end = anti_repeat.anti_repeat_end
            anti_repeat_direction = anti_repeat.anti_repeat_strand

            if anti_repeat_direction == "reversed":
                anti_repeat_upstream_coordinates = (hit_end, hit_end+5)
                anti_repeat_upstream = self.reverse_com(self.dna[hit_end-1:hit_end+4])

                downstream_region_coordinates = (hit_start - 150, hit_start)
                downstream_region = self.reverse_com(self.dna[(hit_start-151):hit_start-1])

                signal_coordinates = self.get_subseqs_by_signal(downstream_region, ("A", "T"))

            else:
                anti_repeat_upstream_coordinates = (hit_start-5, hit_start)
                anti_repeat_upstream = self.dna[(hit_start-6):(hit_start-1)]

                downstream_region_coordinates = (hit_end, hit_end+150)
                downstream_region = self.dna[hit_end - 1:hit_end + 149]

                signal_coordinates = self.get_subseqs_by_signal(downstream_region, ("A", "T"))

            complete_tracr_seq = anti_repeat_seq + downstream_region
            tracr_rna = AntiRepeatTracrRNA(anti_repeat,
                                           complete_tracr_seq,
                                           anti_repeat_upstream,
                                           anti_repeat_seq,
                                           downstream_region,
                                           anti_repeat_upstream_coordinates,
                                           hit_start,
                                           hit_end,
                                           downstream_region_coordinates,
                                           signal_coordinates,
                                           anti_repeat_direction)

            self.list_tracr_rna.append(tracr_rna)

    def output(self):
        return self.list_tracr_rna
