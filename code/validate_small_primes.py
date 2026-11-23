"""Validate the description of rational modular units on X_1(p) against known values.

Run with Sage's python from the repository root:
    python3 code/validate_small_primes.py [largest prime] [exponent bound]
For each prime 11 <= p <= bound we minimize the degree of a nonconstant unit
prod_a H_a^{e_a} (Streng, Theorem 1.2) with |e_a| <= B, using the order formula at the
cusps, and compare the minimum with the Q-gonality of X_1(p) for p <= 40 (Derickx and
van Hoeij, where it is attained by a modular unit) and with their upper bound for p > 40.
This is a floating-point MILP (HiGHS via scipy), used only as a sanity check of the model:
it is a test on known positives, not part of the proof.
"""
import sys
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from sage.all import prime_range

# Derickx and van Hoeij, https://www.math.fsu.edu/~hoeij/files/X1N/gonality
# (exact Q-gonality for N <= 40, upper bound by an explicit modular unit for N > 40).
dvh = {11: 2, 13: 2, 17: 4, 19: 5, 23: 7, 29: 11, 31: 12, 37: 18, 41: 22, 43: 24,
       47: 29, 53: 37, 59: 46, 61: 49, 67: 58, 71: 66, 73: 70, 79: 82, 83: 90,
       89: 104, 97: 123, 101: 133, 103: 139, 107: 150, 109: 156, 113: 167}


def minimal_unit_degree(p, B):
    h = (p - 1) // 2
    # Orders times 12p: 12p * (p/2) B_2(r/p) = 6 r (r - p) + p^2.
    W = np.array([[6 * (a * c % p) * ((a * c % p) - p) + p * p for a in range(1, h + 1)]
                  for c in range(1, h + 1)], dtype=float)
    # variables: e_1..e_h, s, T, z_1..z_h, w
    n = 2 * h + 3
    ie, i_s, i_t, iz, iw = 0, h, h + 1, h + 2, 2 * h + 2
    cost = np.zeros(n)
    cost[iz:iz + h] = 1
    cost[iw] = h
    rows, lo, hi = [], [], []

    def add(row, l, u):
        rows.append(row); lo.append(l); hi.append(u)

    r = np.zeros(n); r[ie:ie + h] = 1; r[i_s] = -12; add(r, 0, 0)          # sum e = 12 s
    r = np.zeros(n); r[ie:ie + h] = [a * a % p for a in range(1, h + 1)]
    r[i_t] = -p; add(r, 0, 0)                                            # square sum = p T
    for c in range(h):                                                   # 12p z_c >= 12p v_c
        r = np.zeros(n); r[ie:ie + h] = W[c]; r[iz + c] = -12 * p; add(r, -np.inf, 0)
    r = np.zeros(n); r[iw] = 1; r[i_s] = -1; add(r, 0, np.inf)           # w >= s
    lb = np.full(n, -np.inf); ub = np.full(n, np.inf)
    lb[ie:ie + h] = -B; ub[ie:ie + h] = B
    lb[ie] = 1                                     # a diamond and inversion give e_1 >= 1
    lb[iz:iz + h] = 0; lb[iw] = 0
    integrality = np.zeros(n); integrality[ie:ie + h + 2] = 1
    res = milp(cost, constraints=LinearConstraint(np.array(rows), lo, hi),
               integrality=integrality, bounds=Bounds(lb, ub),
               options={"time_limit": 600, "mip_rel_gap": 0})
    e = np.round(res.x[ie:ie + h]).astype(int)
    return round(res.fun), res.status, int(np.abs(e).max())


if __name__ == "__main__":
    top = int(sys.argv[1]) if len(sys.argv) > 1 else 83
    B = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    bad = 0
    for p in prime_range(11, top + 1):
        degree, status, largest = minimal_unit_degree(p, B)
        ref = dvh[p]
        if p <= 40:
            ok = degree == ref
            note = f"gon_Q X_1({p}) = {ref}"
        else:
            ok = degree <= ref
            note = f"Derickx-van Hoeij unit of degree {ref}"
        flag = "" if status == 0 else " (time limit, not proven optimal)"
        bad += not ok
        print(f"p = {p:3d}: minimal unit degree {degree:4d}{flag}, max |e_a| = {largest}; {note}"
              f" {'OK' if ok else 'MISMATCH'}", flush=True)
    print("VALIDATION PASSED" if bad == 0 else f"{bad} MISMATCHES")
