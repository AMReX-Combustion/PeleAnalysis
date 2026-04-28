#!/usr/bin/env python3
"""
check_plt_range.py <plt_dir> <var_name> <expected_min> <expected_max> [<tol>]

Reads every cell of <var_name> across all AMR levels and checks that the
observed minimum ≈ expected_min and observed maximum ≈ expected_max within
<tol> (default: 1e-10).  Exits 0 on success, 1 on failure.
"""
import sys
import os
import re
import struct


def parse_fab_hdr(line):
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
    if len(sys.argv) < 5:
        sys.exit(
            "Usage: check_plt_range.py <plt_dir> <var_name>"
            " <expected_min> <expected_max> [tol]"
        )

    plt_dir = sys.argv[1]
    var_name = sys.argv[2]
    expected_min = float(sys.argv[3])
    expected_max = float(sys.argv[4])
    tol = float(sys.argv[5]) if len(sys.argv) > 5 else 1e-10

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

    obs_min = float("inf")
    obs_max = float("-inf")
    total_cells = 0

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
            obs_min = min(obs_min, min(vals))
            obs_max = max(obs_max, max(vals))
            total_cells += n
        lev += 1

    if total_cells == 0:
        sys.exit(f"No cells found in {plt_dir}")

    ok = True
    errs = []
    if abs(obs_min - expected_min) > tol:
        ok = False
        errs.append(
            f"min mismatch: got {obs_min:.6g}, expected {expected_min:.6g}"
            f" (err={abs(obs_min-expected_min):.3e} > tol={tol})"
        )
    if abs(obs_max - expected_max) > tol:
        ok = False
        errs.append(
            f"max mismatch: got {obs_max:.6g}, expected {expected_max:.6g}"
            f" (err={abs(obs_max-expected_max):.3e} > tol={tol})"
        )

    summary = (
        f"Checked {total_cells} cells: min={obs_min:.6g}, max={obs_max:.6g}"
        f"  (expected [{expected_min}, {expected_max}] ±{tol})"
    )
    if errs:
        print(summary + "  FAIL: " + "; ".join(errs))
    else:
        print(summary)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
