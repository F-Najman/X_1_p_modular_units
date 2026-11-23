"""Exact necessary conditions on a degree <= 380 modular unit on X_1(191).

Run with sage -python cusp_bounds.py. In the notation of Streng, Theorem 1.2,
e is the vector of integral Siegel exponents, v=Ae is the vector of orders at
the rational cusps, and s=sum(e)/12 is the order at the degree-95 cusp.
"""
from sage.all import QQ, matrix, primitive_root, vector


p, h, d = 191, 95, 380


def fold(x):
    return min(x % p, (-x) % p)


def cusp_order(c, a):
    r = c * a % p
    return QQ(6 * r * r - 6 * r * p + p * p) / (12 * p)


A = matrix(QQ, h, h, lambda i, j: cusp_order(i + 1, j + 1))
unit_vector = vector(QQ, [1] + [0] * (h - 1))
column = A.solve_right(unit_vector)
M = matrix(QQ, h, h, lambda i, j: column[fold((i + 1) * (j + 1)) - 1])
assert A * M == matrix.identity(QQ, h)
assert set(sum(row) for row in A.rows()) == {-QQ(h) / 12}


def upper(s, row):
    # Invert the unit first to arrange s >= 0. The positive mass of v is at
    # most d-h*s and the negative mass at most d. The two measures may be
    # padded by equal amounts at the same coordinate to attain these masses.
    return (d - h * s) * max(row) - d * min(row)


for s in range(d // h + 1):
    low = -upper(s, -M[0]).floor()
    high = upper(s, M[0]).floor()
    print("cusp order", s, "individual exponents", low, high)
    assert (low, high) == [(-2, 2), (-1, 1), (-1, 1), (0, 1), (0, 1)][s]

# If e_1=2 and s=0, the two inequalities force -1 < e_2 < 0.
# Diamonds and inversion cover every exponent of absolute value two.
assert upper(0, 16*M[0] - M[1]) < 33
assert upper(0, 4*M[0] + M[1]) < 8

def top_sum(values, length):
    return sum(sorted(values, reverse=True)[:length])


# At s=4 the nonzero exponents form a 48-element set S. All rational cusp
# orders are nonpositive, and their negatives n_c have total d. Put
# L_c = -sum_{i in S} M[i,c]. Then 48=sum_c n_c L_c.
k = 48
peak = top_sum(-column, k)
average = QQ(k) / d
pair = max(top_sum(-column - M.column(j), k) for j in range(1, h))
second = pair - average
assert 0 < second < average < peak
# Move a cusp at which L_c is largest to c=1 by a diamond. At other cusps
# L_c <= second. The identity above forces a pole at c=1 of order >= 339.
needed = (k - d * second) / (peak - second)
allowed = top_sum(-A.column(0), k)
assert needed.ceil() == 339 and allowed.floor() == 285
print("s=4, one pole would have order at least", needed.ceil(),
      "and at most", allowed.floor())

# Sums of arbitrary exponent coordinates give further necessary conditions.
for s in range(3):
    negative_counts = []
    for n in range(h + 1):
        m = n + 12*s
        if m + n > h:
            break
        upper_pos = (d - h*s)*top_sum(column, m) + d*top_sum(-column, m)
        upper_neg = (d - h*s)*top_sum(-column, n) + d*top_sum(column, n)
        if m <= upper_pos and n <= upper_neg:
            negative_counts.append(n)
    print("s =", s, "negative exponent count <=", max(negative_counts))

# A nontrivial stabilizer contains a subgroup of order 5 or 19. Invariance
# of the divisor makes e constant on its orbits, since A is invertible.
# The orbit size divides 12*s and hence s=0. All cusp orders repeat in
# equal-sized blocks. Invert the resulting smaller valuation matrix.
g = int(primitive_root(p))
for size in (5, 19):
    count = h // size
    orbits = [[fold(pow(g, i + count*j, p)) for j in range(size)]
              for i in range(count)]
    assert set().union(*map(set, orbits)) == set(range(1, h+1))
    assert all(len(set(orbit)) == size for orbit in orbits)
    folded = matrix(QQ, count, count,
                    lambda i, j: sum(cusp_order(orbits[i][0], a) for a in orbits[j]))
    inverse = folded.inverse()
    cap = max((max(row) - min(row)) * QQ(d)/size for row in inverse.rows())
    assert cap < 1  # Every exponent would have to be zero.
    print("diamond subgroup", size, "exponent bound below", cap.n(digits=7))
