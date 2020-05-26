def diff(filepath_a, filepath_b):
    """
    Difference between two files.
    """
    from difflib import ndiff

    with open(filepath_a, "r") as file_a:
        a = file_a.readlines()

    with open(filepath_b, "r") as file_b:
        b = file_b.readlines()

    return "".join([line.expandtabs(1) for line in ndiff(a, b)])
