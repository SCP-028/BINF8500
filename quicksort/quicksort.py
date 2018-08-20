import sys


def QuickSort(lol, start, stop):
    """Implements the QuickSort algorithm.

    Parameters
    ----------
        lol  : list[list[str]]
        start: int
        stop : int

    Returns
    -------
        A `lol` sorted alphabetically by the first item in each list.
    """
    if start < stop:
        i = start  # i is going to be the correct position for pivot
        pivot = lol[stop][1]
        for j in range(start, stop):
            if lol[j][1] <= pivot:
                lol[i], lol[j] = lol[j], lol[i]
                i += 1
        lol[i], lol[stop] = lol[stop], lol[i]  # put pivot at the right position

        if i > len(lol) / 2:
            lol = QuickSort(lol, i + 1, stop)
            lol = QuickSort(lol, start, i - 1)
        else:
            lol = QuickSort(lol, start, i - 1)
            lol = QuickSort(lol, i + 1, stop)

    return lol


if __name__ == "__main__":
    args = sys.argv[1:]
    print(f"From file {args[0]} to file {args[1]}...")
    with open(args[0], "r") as f:
        fastq = f.readlines()
    fastq = [fastq[i:i + 4] for i in range(0, len(fastq), 4)]

    fastq = QuickSort(fastq, 0, len(fastq) - 1)
    with open(args[1], "w") as f:
        f.write("".join([x for item in fastq for x in item]))
