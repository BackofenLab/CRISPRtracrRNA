def candidate_ranking(candidate_csv_file, separator, dict_weights):
    with open(candidate_csv_file, 'r') as csv_file:
        lines = csv_file.readlines()
        header, info_lines = lines[0], lines[1:]
        header_values = header.strip().split(separator)
        chunks =


def split_into_ch