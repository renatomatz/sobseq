import re
import argparse
import numpy as np


ROWS = 21200
MAX_COLS = 21


def clean(line):
    return re.sub("\s+", ',', line.replace('\n', '')).split(',')


def to_arr(cl):
    return np.array(cl + ['0']*(MAX_COLS-len(cl)), dtype=np.int32)


def rw(name_in, name_out=None, n_rows=ROWS, max_cols=MAX_COLS):
    if name_out is None:
        name_out = name_in
    arr = np.empty([n_rows, max_cols])
    with open(name_in, 'r') as f:
        for i, line in enumerate(f.readlines()):
            arr[i, :] = to_arr(clean(line))
    np.savetxt(name_out, arr, delimiter=', ', fmt="%d")


parser = argparse.ArgumentParser(
    description="Process irregular files into a format "
                "easier for Fortran to read"
)
parser.add_argument("file_names", metavar="names", type=str, nargs='+',
                    help="Names of the files to be processed")
parser.add_argument("--rows", type=int, default=ROWS,
                    help=f"Number of rows being read, defaults to {ROWS}")
parser.add_argument("--max_cols", type=int, default=MAX_COLS,
                    help=f"Maximum number of columns being read, defaults to "
                         f"{MAX_COLS}")


if __name__ == "__main__":
    # TODO: make complete argument parsing

    ns = parser.parse_args()

    for i, name in enumerate(ns.file_names):
        rw(name, f"clean_{name}",
           n_rows=ns.rows,
           max_cols=ns.max_cols)
