"""Cross-check of the 56 integer programs with a second solver (HiGHS, floating point).

Usage (from the repository root):  python3 code/highs_crosscheck.py logs/checks/highs.tsv [workers]
The constraints are rebuilt here from the description in the paper, without using
unit_model.py or the stored MPS files, for every case of checks.tsv. HiGHS works in
floating-point arithmetic, so this is a consistency check, not part of the proof.
"""
import sys
import time
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from scipy.optimize import milp, LinearConstraint, Bounds

p, h, D = 191, 95, 380
root = Path(__file__).resolve().parent.parent
cycle = [min(pow(2, i, p), p - pow(2, i, p)) for i in range(h)]
# 12p times the order of H_a at the cusp <c>infinity: 6 r (r - p) + p^2 with r = ac mod p.
W = np.array([[6 * (a * c % p) * ((a * c % p) - p) + p * p for a in range(1, h + 1)]
              for c in range(1, h + 1)], dtype=float)


def feasible(s, gap, endpoint, fixed):
    n = 2 * h + 1                       # e_1..e_h, z_1..z_h, T
    E, Z, T = 0, h, 2 * h
    rows, lo, hi = [], [], []

    def add(coeffs, l, u):
        r = np.zeros(n)
        for j, x in coeffs:
            r[j] += x
        rows.append(r); lo.append(l); hi.append(u)

    add([(E + a, 1) for a in range(h)], 12 * s, 12 * s)
    add([(E + a, (a + 1) ** 2 % p) for a in range(h)] + [(T, -p)], 0, 0)
    for c in range(h):
        add([(E + a, W[c, a]) for a in range(h)] + [(Z + c, -12 * p)], -np.inf, 0)
    add([(Z + c, 1) for c in range(h)], -np.inf, D - h * s)
    lb = np.zeros(n); ub = np.zeros(n)
    lb[E:E + h] = -1 if s <= 2 else 0; ub[E:E + h] = 1
    lb[E] = ub[E] = 1                   # e_1 = 1
    lb[Z:Z + h] = 0; ub[Z:Z + h] = D - h * s
    lb[T], ub[T] = -h, h
    for a, x in fixed:
        lb[E + a - 1] = ub[E + a - 1] = x
    if gap is not None:
        b = cycle[gap]
        value = 1 if endpoint is None else endpoint
        lb[E + b - 1] = ub[E + b - 1] = value
        for t in range(1, gap):
            for i in range(h):
                x, y = cycle[i] - 1, cycle[(i + t) % h] - 1
                if endpoint is None:    # two positive exponents at distance t < gap
                    add([(E + x, 1), (E + y, 1)], -np.inf, 1)
                else:                   # two nonzero exponents at distance t < gap
                    for sx, sy in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
                        add([(E + x, sx), (E + y, sy)], -np.inf, 1)
    integrality = np.ones(n); integrality[Z:Z + h] = 0
    start = time.time()
    res = milp(np.zeros(n), constraints=LinearConstraint(np.array(rows), lo, hi),
               integrality=integrality, bounds=Bounds(lb, ub),
               options={"time_limit": 3600})
    return res.status, res.message, time.time() - start


def run(row):
    case, order, gap, endpoint, fixed = row[:5]
    fixed = [tuple(map(int, t.split("="))) for t in fixed.split(",")] if fixed else []
    status, message, seconds = feasible(int(order), int(gap) if gap else None,
                                        int(endpoint) if endpoint else None, fixed)
    verdict = "infeasible" if status == 2 else f"status {status}: {message}"
    print(case, verdict, f"{seconds:.1f}s", flush=True)
    return case, verdict, seconds


if __name__ == "__main__":
    out = Path(sys.argv[1])
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    rows = [line.split("\t") for line in (root / "checks.tsv").read_text().splitlines()[1:]]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(run, rows))
    out.write_text("".join(f"{c}\t{v}\t{t:.1f}\n" for c, v, t in results))
    bad = [c for c, v, _ in results if v != "infeasible"]
    print("HIGHS: ALL 56 CASES INFEASIBLE" if not bad else f"HIGHS: NOT INFEASIBLE: {bad}")
