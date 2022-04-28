HEADER_ANTI_REPEAT = 'intarna_anti_repeat_interaction_interval'
HEADER_TAIL_MODEL_HIT_LOCATION = 'tail_model_hit_location'
HEADER_BEST_TERMINATOR_LOCATION = 'best_terminator_location'
HEADER_CRISPR_ORIENTATION = 'crispr_array_orientation'


def consistency_score_maker(csv_file, separator, output_file):
    with open(csv_file, 'r') as f:
        lines = f.readlines()
        header, info_lines = lines[0], lines[1:]

        new_header = header.strip() + f"{separator}" + f"{separator}".join(["consistency_anti_repeat_tail", "consistency_tail_terminator", "consistency_terminator_orientation"]) + '\n'
        lines_new = []

        index_anti_repeat = get_header_column_index(header, separator, HEADER_ANTI_REPEAT)
        index_tail = get_header_column_index(header, separator, HEADER_TAIL_MODEL_HIT_LOCATION)
        index_terminator = get_header_column_index(header, separator, HEADER_BEST_TERMINATOR_LOCATION)
        index_orientation = get_header_column_index(header, separator, HEADER_CRISPR_ORIENTATION)

        for line in info_lines:
            interval_anti_repeat = get_interval_anti_repeat(index_anti_repeat, line, separator)
            interval_tail = get_interval_tail(index_tail, line, separator)
            interval_terminator = get_interval_terminator(index_terminator, line, separator)
            orientation = get_orientation(index_orientation, line, separator)
            consistency_anti_repeat_tail = get_consistency_anti_repeat(interval_anti_repeat, interval_tail)
            consistency_tail_terminator = get_consistency_tail_terminator(interval_tail, interval_terminator)
            consistency_terminator_orientation = get_consistency_terminator_orientation(orientation, interval_terminator)
            new_line = line.strip() + separator + str(consistency_anti_repeat_tail) + separator + str(consistency_tail_terminator) + separator + str(consistency_terminator_orientation) + '\n'
            lines_new.append(new_line)

    with open(output_file, 'w') as f:
        f.write(new_header)
        f.writelines(lines_new)


def get_header_column_index(header, separator, header_name):
    header = header.strip().split(separator)
    for i in range(len(header)):
        if header[i] == header_name:
            return i


def get_orientation(orientation_index, line, separator):
    line_elements = line.strip().split(separator)
    crispr_orientation = line_elements[orientation_index]
    return crispr_orientation


def get_interval_anti_repeat(anti_repeat_index, line, separator):
    line_elements = line.strip().split(separator)
    anti_repeat_interval = line_elements[anti_repeat_index]
    if anti_repeat_interval == 'NA':
        return None
    if anti_repeat_interval in ("", " "):
        return None
    else:
        start = int(anti_repeat_interval.split('_')[0])
        end = int(anti_repeat_interval.split('_')[1])
        return start, end


def get_interval_tail(tail_index, line, separator):
    line_elements = line.strip().split(separator)
    tail_model_hit_location = line_elements[tail_index]
    if tail_model_hit_location == 'NA':
        return None
    else:
        start = int(tail_model_hit_location.split('_')[0])
        end = int(tail_model_hit_location.split('_')[1])
        return start, end


def get_interval_terminator(terminator_index, line, separator):
    line_elements = line.strip().split(separator)
    best_terminator_location = line_elements[terminator_index]
    if best_terminator_location == 'NA':
        return None
    else:
        best_terminator_location=best_terminator_location.replace("(", "").replace(")", "").replace("- ", "_")
        start = int(best_terminator_location.split('_')[0])
        end = int(best_terminator_location.split('_')[1])
        return start, end


def overlap(interval1, interval2, overlap_size):
    if interval1 is None or interval2 is None:
        return False
    else:
        start1, end1 = interval1
        start2, end2 = interval2
        if start1 <= start2 <= end1:
            if end1 - start2 >= overlap_size:
                return True
            else:
                return False
        if start2 <= end1 <= end2:
            if end2 - start1 >= overlap_size:
                return True
            else:
                return False


def get_consistency_anti_repeat(interval_anti_repeat, interval_tail):
    if overlap(interval_anti_repeat, interval_tail, 10):
        return 0
    else:
        return 1


def get_consistency_tail_terminator(interval_tail, interval_terminator):
    if overlap(interval_tail, interval_terminator, 10):
        return 0
    else:
        return 1


def get_consistency_terminator_orientation(orientation, interval_terminator):
    if orientation == 'Forward':
        if interval_terminator:
            return 1
        else:
            return 0
    else:
        if interval_terminator:
            return 0
        else:
            return 1

