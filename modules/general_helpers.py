def reverse_com(seq):
    dict_complementary = {"A": "T", "G": "C", "T": "A", "C": "G", "N": "N"}
    new_seq = "".join([dict_complementary.get(x, x) for x in seq])
    new_seq = new_seq[::-1]
    return new_seq