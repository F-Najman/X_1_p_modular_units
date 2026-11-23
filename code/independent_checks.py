"""Independent checks of the facts about modular units on X_1(191) used in the paper.

Run with Sage's python from the repository root:  python3 code/independent_checks.py
This file does not import any of the other scripts. It recomputes, in exact rational
arithmetic:
  Claim 1: the order matrix A is invertible, its columns sum to -h/12, and the unit
           H_4 H_7 / (H_1 H_8) has degree 478;
  Claim 2: the bounds on the exponents e_a of a unit of degree <= 380 (Proposition 5.2);
  Claim 3: an exponent equal to 2 is impossible when s = 0 (Proposition 5.2);
  Claim 4: the multiplication-by-2 cycle on (Z/191Z)^*/{+-1} has length 95.
"""
from itertools import product
from sage.all import QQ, ZZ, matrix, vector, floor

p, h, D = 191, 95, 380


def B2(x):
    return x * x - x + QQ(1) / 6


def frac(x):
    return x - floor(x)


def order_matrix(p):
    h = (p - 1) // 2
    # Row c, column a: the order of H_a at the cusp <c>infinity, in the local
    # parameter q^(1/p), namely (p/2) B_2({ac/p}).
    return matrix(QQ, h, h, lambda c, a: QQ(p) / 2 * B2(frac(QQ((a + 1) * (c + 1)) / p)))


A = order_matrix(p)

# Claim 1
assert A.det() != 0
assert all(sum(A.column(a)) == -QQ(h) / 12 for a in range(h))
e = vector(QQ, h)
for a, x in [(4, 1), (7, 1), (1, -1), (8, -1)]:
    e[a - 1] = x
assert sum(e) == 0 and sum((a + 1) ** 2 * e[a] for a in range(h)) % p == 0
v = A * e
assert all(x in ZZ for x in v) and sum(v) == 0
assert sum(max(x, 0) for x in v) == 478
print("Claim 1 successfully verified: A invertible, column sums -95/12, deg H4H7/(H1H8) = 478")

# Claim 2. For s >= 0 and a unit of degree <= D, the vector v of orders at the rational
# cusps lies in Q_s = {v : sum max(v_c,0) <= D - h s, sum v_c = -h s}, and e = M v.
# The maximum of a linear form r.v over Q_s is (D - h s) max(r) - D min(r).
M = A.inverse()


def lp_max(r, s):
    return (D - h * s) * max(r) - D * min(r)


expected = {0: (-2, 2), 1: (-1, 1), 2: (-1, 1), 3: (0, 1), 4: (0, 1)}
for s in range(5):
    lows, highs = set(), set()
    for a in range(h):
        r = M.row(a)
        highs.add(floor(lp_max(r, s)))
        lows.add(-floor(lp_max(-r, s)))
    assert (min(lows), max(highs)) == expected[s], (s, lows, highs)
    print(f"  s = {s}: every exponent lies in [{min(lows)}, {max(highs)}]")
assert D - h * 5 < 0  # s >= 5 is impossible, as the degree is at least h s
print("Claim 2 successfully verified")

# Claim 3. Suppose s = 0 and e_1 = 2. We compute the exact maximum and minimum of e_2
# over {v in Q_0 : (M v)_1 = 2} by LP duality: for a linear form r and the constraint
# m.v = 2, the maximum of r.v equals min over lambda of lp_max(r - lambda m, 0) + 2 lambda,
# a convex piecewise linear function of lambda whose minimum is attained at a breakpoint.


def constrained_max(r, m):
    breakpoints = {(r[i] - r[j]) / (m[i] - m[j])
                   for i in range(h) for j in range(h) if m[i] != m[j]}
    return min(lp_max(r - lam * m, 0) + 2 * lam for lam in breakpoints)


m1, m2 = M.row(0), M.row(1)
top = constrained_max(m2, m1)
bottom = -constrained_max(-m2, m1)
assert -1 < bottom and top < 0, (bottom, top)
print(f"Claim 3 successfully verified: if s = 0 and e_1 = 2 then {bottom.n(digits=6)} <= e_2 <= {top.n(digits=6)}")

# Claim 4
cycle, x = [], 1
for _ in range(h):
    cycle.append(min(x, p - x))
    x = 2 * x % p
assert sorted(cycle) == list(range(1, h + 1)) and min(x, p - x) == 1
print("Claim 4 successfully verified: 2 generates (Z/191Z)^*/{+-1}")
