#!/usr/bin/env python3
"""
check_plt_coord.py <plt_dir> <expected_coord_sys>

Reads the AMReX plotfile Header and asserts that the coord_sys field equals
<expected_coord_sys> (0 = Cartesian, 1 = cylindrical/RZ).

The coord_sys integer appears in the Header after the per-level cell-size
lines.  Formula: coord_line_idx = 11 + ncomp + finest_level (0-indexed).

Exits 0 on success, 1 on mismatch or missing file.
"""
import sys
import os


def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: check_plt_coord.py <plt_dir> <expected_coord_sys>")

    plt_dir = sys.argv[1]
    expected = int(sys.argv[2])

    hdr_path = os.path.join(plt_dir, "Header")
    if not os.path.isfile(hdr_path):
        print(f"Header not found: {hdr_path}")
        sys.exit(1)

    with open(hdr_path) as f:
        lines = f.readlines()

    ncomp = int(lines[1].strip())
    finest_level = int(lines[4 + ncomp].strip())
    coord_line_idx = 11 + ncomp + finest_level

    if coord_line_idx >= len(lines):
        print(f"Header too short (need line {coord_line_idx}, have {len(lines)})")
        sys.exit(1)

    coord_sys = int(lines[coord_line_idx].strip())

    if coord_sys == expected:
        print(f"coord_sys = {coord_sys} (as expected)")
        sys.exit(0)
    else:
        print(f"coord_sys = {coord_sys}, expected {expected}")
        sys.exit(1)


if __name__ == "__main__":
    main()
