import re


def index_header_names(header_values):
    dict_indexes = {}
    for i in range(len(header_values)):
        dict_indexes[header_values[i]] = i
    return dict_indexes


def create_regex_pattern_from_dna_mask(mask_sequence):
    dict_extended = {
        'A': ["A"],
        'C': ["C"],
        'G': ["G"],
        'T': ["T"],
        'R': ["A", "G"],
        'Y': ["C", "T"],
        'M': ["A", "C"],
        'K': ["G", "T"],
        'S': ["G", "C"],
        'W': ["A", "T"],
        'H': ["A", "C", "T"],
        'B': ["C", "G", "T"],
        'V': ["A", "C", "G"],
        'D': ["A", "G", "T"],
        'N': ["A", "C", "G", "T"]
    }
    regex_pattern = ""
    for letter in mask_sequence:
        regex_pattern += "[" + "|".join(dict_extended[letter]) + "]"

    return regex_pattern


def search_for_pattern_with_re(mask_sequence, text):
    pattern = create_regex_pattern_from_dna_mask(mask_sequence)
    result = re.search(pattern, text)
    if result:
        interval = result.span()
        interval = interval[0] + 1, interval[1] + 1
        match = result.group()
        return str(interval)[1:-1].replace(", ", "-"), match
    else:
        return "NA", "NA"


def search_two_patterns_with_interval(pattern1, pattern2, low_boundary, high_boundary, text):
    starts_first = [m.start() for m in re.finditer(pattern1, text)]
    starts_second = [m.start() for m in re.finditer(pattern2, text)]
    for start_first in starts_first:
        for start_second in starts_second:
            distance = start_second - start_first + len(pattern1)
            if low_boundary <= distance <= high_boundary:
                interval_first = str(start_first+1) + "-" + str(start_first + len(pattern1)+1)
                interval_second = str(start_second+1) + "-" + str(start_second + len(pattern1)+1)
                return interval_first, interval_second

    return "NA", "NA"


def main():
    mask_sequence = "CCYCCNNNNNNGGRGG"
    text = "AAAAAAAAAAAAACCCCCCAAAAAAGGAGGCCCCCAAAACCCCCAAGGAGGAAAAAAAAAAAAACCCCCCAAAAAAGGAGGCCCCCAAAACCCCCAAGGAGG"
    result = search_for_pattern_with_re(mask_sequence, text)
    print(result)

    mask_sequence = "CCYCCNNNNNNGGRGG"
    text = "N"
    result = search_for_pattern_with_re(mask_sequence, text)
    print(result)

    pattern1 = "TTT"
    pattern2 = "CTTTC"
    low_boundary = 12
    high_boundary = 28
    text = "TTTAAAAAAAAAAAAAAAAAAAAAACTTTC"
    result = search_two_patterns_with_interval(pattern1, pattern2, low_boundary, high_boundary, text)
    print(result)


def anti_repeat_search_type_v(csv_file_name_old, delimiter, csv_file_name_new):
    with open(csv_file_name_old, 'r') as csv_file:
        lines = csv_file.readlines()
        header, info_lines = lines[0], lines[1:]
        header_values = header.strip().split(delimiter)
        dict_indexes = index_header_names(header_values)
        all_hits_anti_repeat1 = []
        all_hits_anti_repeat2 = []
        for line in info_lines:
            tracr_hit_sequence = line.split(delimiter)[dict_indexes["hit_sequence"]]
            hits_anti_repeat_1 = search_two_patterns_with_interval("TTT", "CTTTC", 12, 28, tracr_hit_sequence)
            hits_anti_repeat_2 = search_for_pattern_with_re("CCYCCNNNNNNGGRGG", tracr_hit_sequence)
            all_hits_anti_repeat1.append(hits_anti_repeat_1)
            all_hits_anti_repeat2.append(hits_anti_repeat_2)

    with open(csv_file_name_new, 'w') as csv_file:
        new_header = header.strip() + delimiter + "region_anti_repeat1 TTT" + delimiter + "region_anti_repeat1 CTTTC" + delimiter + "interval_motif_before_anti_repeat2" + delimiter + "sequence_motif_before_anti_repeat2" "\n"
        csv_file.write(new_header)

        for line, info_anti_repeat1, info_anti_repeat2 in zip(info_lines, all_hits_anti_repeat1, all_hits_anti_repeat2):
            new_line = line.strip() + delimiter + info_anti_repeat1[0] + delimiter + info_anti_repeat1[1] + delimiter + info_anti_repeat2[0] + delimiter + info_anti_repeat2[1] + "\n"
            csv_file.write(new_line)


if __name__ == "__main__":
    main()
