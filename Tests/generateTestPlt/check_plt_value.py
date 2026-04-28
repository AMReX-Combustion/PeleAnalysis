#!/usr/bin/env python3
"""
check_plt_value.py <plt_dir> <var_name> <expected_value> [<tolerance>]

Reads every cell of variable <var_name> across all AMR levels of an AMReX
plotfile and asserts that each value equals <expected_value> within <tolerance>
(default: 1e-10).  Exits 0 on success, 1 on failure, 2 on usage error.

The plotfile must be in the standard AMReX binary FArrayBox format produced by
WriteMultiLevelPlotfile.  Component-major storage is assumed (all values of
component k before component k+1).
"""
import sys
import os
import re
import struct
import glob


def parse_fab_hdr(line):
    """Return (ncomp, lo, hi, little_endian) from an AMReX FAB header line.

    AMReX encodes endianness as (8, (8 7 6 5 4 3 2 1)) for little-endian and
    (8, (1 2 3 4 5 6 7 8)) for big-endian.
    """
    little = "8 7 6 5 4 3 2 1" in line
    triples = re.findall(r"\((-?\d+),(-?\d+),(-?\d+)\)", line)
    if triples:
        lo = tuple(int(x) for x in triples[0])
        hi = tuple(int(x) for x in triples[1])
    else:
        pairs = re.findall(r"\((-?\d+),(-?\d+)\)", line)
        if len(pairs) < 2:
            raise ValueError(f"Cannot parse box from FAB header: {line!r}")
        lo = tuple(int(x) for x in pairs[0])
        hi = tuple(int(x) for x in pairs[1])
    ncomp = int(line.strip().split()[-1])
    return ncomp, lo, hi, little


def npts(lo, hi):
    n = 1
    for a, b in zip(lo, hi):
        n *= b - a + 1
    return n


def main():
    if len(sys.argv) < 4:
        sys.exit(
            "Usage: check_plt_value.py <plt_dir> <var_name> <expected> [tol]"
        )

    plt_dir = sys.argv[1]
    var_name = sys.argv[2]
    expected = float(sys.argv[3])
    tol = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-10

    hdr_path = os.path.join(plt_dir, "Header")
    if not os.path.isfile(hdr_path):
        sys.exit(f"Header not found: {hdr_path}")
    with open(hdr_path) as f:
        hdr_lines = f.readlines()

    ncomp = int(hdr_lines[1].strip())
    var_names = [hdr_lines[2 + i].strip() for i in range(ncomp)]
    if var_name not in var_names:
        sys.exit(
            f"Variable '{var_name}' not in plotfile.  Available: {var_names}"
        )
    comp_idx = var_names.index(var_name)

    total_cells = 0
    max_err = 0.0
    lev = 0
    while True:
        lev_dir = os.path.join(plt_dir, f"Level_{lev}")
        if not os.path.isdir(lev_dir):
            break

        cell_h_path = os.path.join(lev_dir, "Cell_H")
        with open(cell_h_path) as f:
            cell_h = f.read()

        for m in re.finditer(r"FabOnDisk:\s+(\S+)\s+(\d+)", cell_h):
            fab_file = os.path.join(lev_dir, m.group(1))
            byte_offset = int(m.group(2))
            with open(fab_file, "rb") as f:
                f.seek(byte_offset)
                fab_hdr = f.readline().decode("latin-1")
                _, lo, hi, little = parse_fab_hdr(fab_hdr)
                n = npts(lo, hi)
                endian = "<" if little else ">"
                f.seek(comp_idx * n * 8, 1)
                raw = f.read(n * 8)
                if len(raw) < n * 8:
                    sys.exit(
                        f"Truncated FAB data in {fab_file} at offset {byte_offset}"
                    )
                vals = struct.unpack(endian + "d" * n, raw)
            err = max(abs(v - expected) for v in vals)
            max_err = max(max_err, err)
            total_cells += n
        lev += 1

    if total_cells == 0:
        sys.exit(f"No cells found in {plt_dir}")

    print(
        f"Checked {total_cells} cells, var='{var_name}', "
        f"expected={expected}, max_err={max_err:.3e}"
    )
    sys.exit(0 if max_err <= tol else 1)


if __name__ == "__main__":
    main()
