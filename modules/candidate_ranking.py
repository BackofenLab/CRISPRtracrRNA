


def candidate_by_acc_number_ranking(candidate_csv_file, separator, dict_weights):
    with open(candidate_csv_file, 'r') as csv_file:
        lines = csv_file.readlines()
        header, info_lines = lines[0], lines[1:]
        header_values = header.strip().split(separator)
        chunks = split_into_chunks(info_lines, separator)


def candidate_ranking(candidate_csv_file, separator, dict_weights, candidate_ranking_file):
    with open(candidate_csv_file, 'r') as csv_file:
        lines = csv_file.readlines()
        header, info_lines = lines[0], lines[1:]
        header_values = header.strip().split(separator)
        dict_indexes = get_column_indexes(header_values, dict_weights)
        pairs_line_score = []
        for line in info_lines:
            pairs_line_score.append((line, score_line(line, dict_weights, dict_indexes, separator)))
        pairs_line_score.sort(key=lambda x: x[1], reverse=True)

    with open(candidate_ranking_file, 'w') as csv_file:
        new_header = header.strip() + separator + 'score'
        csv_file.write(header)
        for line, score in pairs_line_score:
            csv_file.write(line)



def score_line(line, dict_weights, dict_indexes, separator):
    line_elements = line.strip().split(separator)
    score = 0
    for key in dict_weights:
        score += float(line_elements[dict_indexes[key]]) * dict_weights[key]
    return score


def get_column_indexes(header_values, dict_weights):
    dict_indexes = {}
    for key, value in dict_weights.items():
        dict_indexes[key] = header_values.index(value)
    return dict_indexes


def split_into_chunks(lines, separator):
    chunks = []
    current_chunk = []
    current_acc_num = None
    for line in lines:
        values = line.strip().split(separator)
        acc_number = values[0]
        if current_acc_num is None:
            current_acc_num = acc_number
            current_chunk.append(line)
        elif acc_number == current_acc_num:
            current_chunk.append(line)
        else:
            chunks.append(current_chunk)
            current_chunk = [line]
            current_acc_num = acc_number
    if current_chunk:
        chunks.append(current_chunk)
    return chunks
