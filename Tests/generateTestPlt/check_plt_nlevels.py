#!/usr/bin/env python3
"""
check_plt_nlevels.py <plt_dir> <expected_nlevels>

Counts AMR levels (Level_N/ directories containing Cell_H) in <plt_dir>
and compares with expected_nlevels.  Exits 0 on match, 1 on mismatch.
"""
import sys
import os


def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: check_plt_nlevels.py <plt_dir> <expected_nlevels>")

    plt_dir = sys.argv[1]
    expected = int(sys.argv[2])

    if not os.path.isfile(os.path.join(plt_dir, "Header")):
        sys.exit(f"Header not found in {plt_dir}")

    actual = 0
    while os.path.isdir(os.path.join(plt_dir, f"Level_{actual}")):
        actual += 1

    if actual == expected:
        print(f"nlevels = {actual} (as expected)")
        sys.exit(0)
    else:
        print(f"ERROR: expected {expected} levels, found {actual} in {plt_dir}")
        sys.exit(1)


if __name__ == "__main__":
    main()
