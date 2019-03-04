#!/usr/bin/env python


import sys
from math import log


def not_null(l):
    nl = []
    for li in l:
        if li:
            nl += [li]
    return nl


if not len(sys.argv) > 1:
    print("Name of study file required.")
    exit(1)

with open(sys.argv[1]) as f:
    data = f.readlines()
data = [row[:-1] for row in data]

nrows = (len(data) - 2) // 3
if not nrows * 3 + 2 == len(data):
    print("Invalid number of rows in file.")
    exit(1)

out   = []

print(" & ".join([
    r"\( N \)",
    r"\( L_1 \) Error",
    r"\( L_1 \) Order",
    r"\( L_\infty \) Error",
    r"\( L_\infty \) Order"
    ]) + " \\\\\\hline\\hline")
first_row = True
for ir in range(nrows):
    if not first_row:
        prevN   = N
        prevL1Err = L1Err
        prevLInfErr = LInfErr
    N       = int(           data[3 * ir + 0 + 2].split('-') [2])
    L1Err   = float(not_null(data[3 * ir + 1 + 2].split(' '))[1])
    LInfErr = float(not_null(data[3 * ir + 2 + 2].split(' '))[1])
    line = []
    line += ["{}".format(N)]
    line += ["{:.4f}".format(L1Err)]
    line += ["--" if first_row else "{:.2f}".format(
        log(L1Err / prevL1Err) / log(prevN / N)
        )]
    line += ["{:.4f}".format(LInfErr)]
    line += ["--" if first_row else "{:.2f}".format(
        log(LInfErr / prevLInfErr) / log(prevN / N)
        )]
    print(" & ".join(line) + " \\\\\\hline")
    first_row = False




