def reverse_com(seq):
    dict_complementary = {"A": "T", "G": "C", "T": "A", "C": "G", "N": "N"}
    new_seq = "".join([dict_complementary.get(x, x) for x in seq])
    new_seq = new_seq[::-1]
    return new_seq


def string_tuple_to_tuple(string_tuple):
    """
    :param string_tuple:
    :return: tuple
    Example:
    '>>> string_tuple_to_tuple("(1,2)")'
    (1,2)
    """
    return eval(string_tuple)


def distance_between_two_intervals(interval1, interval2):
    """
    :param interval1:
    :param interval2:
    :return: distance between two intervals
    """
    start1, end1 = interval1
    start2, end2 = interval2
    if start1 < start2:
        if end1 < start2:
            return start2 - end1
        else:
            return 0
    else:
        if end2 < start1:
            return start1 - end2
        else:
            return 0


def get_closest_interval(interval, list_intervals):
    """
    :param interval:
    :param list_intervals:
    :return: distance_to_closest, closest interval
    """
    if not interval:
        return "NA", "NA"
    if interval in ("NA", "", " ", "None"):
        return "NA", "NA"
    if not list_intervals:
        return "NA", "NA"
    distance_to_closest = float("inf")
    closest_interval = None
    for interval_ in list_intervals:
        distance = distance_between_two_intervals(interval, interval_)
        if distance < distance_to_closest:
            distance_to_closest = distance
            closest_interval = interval_
    return distance_to_closest, closest_interval
