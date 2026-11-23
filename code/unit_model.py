"""Write an integer-coefficient LP instance for modular units of degree <= 380.

Usage: python unit_model.py 2 > order2.lp
The argument is the order s at the closed cusp of degree 95, after inversion.
To use SCIP 10.1.0 in exact mode, first convert to fixed MPS with
glpsol --lp order2.lp --wmps order2.mps --check. Enable exact solving
before reading the MPS file. The optional gap breaks diamond symmetry.
"""
import argparse


def linear(terms):
    parts = []
    for coefficient, variable in terms:
        if coefficient:
            sign = "-" if coefficient < 0 else "+"
            parts.append(f"{sign} {abs(coefficient)} {variable}")
    return " ".join(parts).removeprefix("+ ")


def write_model(s, gap=None, endpoint=None, reduced=False, fixed=()):
    p, h, d = 191, 95, 380
    indices = range(1, h + 1)
    positive_mass = d - h*s
    e = {a: f"e{a}" for a in indices}
    z = {a: f"z{a}" for a in indices}
    mass = linear((1, z[c]) for c in indices)
    print("Minimize\n cost:", mass)
    print("Subject To")
    print(" mass:", mass, "<=", positive_mass)
    print(" sum_e:", linear((1, e[a]) for a in indices), "=", 12*s)
    def centered(x):
        return (x + h) % p - h
    residues = {a: centered(a*a*pow(2, -1, p)) if reduced else a*a % p for a in indices}
    print(" square_sum:", linear((residues[a], e[a]) for a in indices),
          f"- {p} multiple = 0")
    if gap is not None:
        # Place a shortest gap between successive positive exponents at 1.
        # There are at least 12*s positives if s>0, and at least 6 if s=0
        # by the exhaustive balanced-products test through five factors.
        cycle = [min(pow(2,i,p), (-pow(2,i,p)) % p) for i in range(h)]
        assert len(set(cycle)) == h
        if endpoint is None:
            assert 1 <= gap <= h // (12*s if s else 6)
            right = 1
        else:
            # At s=0 there are at least 12 nonzeros. Inversion normalizes
            # the first sign to +. The shortest nonzero gap is at most 7.
            assert s == 0 and endpoint in (-1, 1) and 1 <= gap <= h//12
            right = endpoint
        print(f" endpoint: {e[cycle[gap]]} = {right}")
        for shift in range(1, gap):
            for i, a in enumerate(cycle):
                b = cycle[(i + shift) % h]
                signs = [(1,1)] if endpoint is None else [(1,1),(1,-1),(-1,1),(-1,-1)]
                for t, (x,y) in enumerate(signs):
                    print(f" gap{shift}_{i}_{t}:", linear([(x,e[a]),(y,e[b])]), "<= 1")
    for c in indices:
        terms = []
        factor = centered(c*c)
        for a in indices:
            r = a*c % p
            coefficient = r*(r-p)//2
            if reduced:
                coefficient -= factor*residues[a]
                assert coefficient % p == 0
                coefficient //= p
            terms.append((coefficient, e[a]))
        if reduced:
            terms.extend([(factor,"multiple"),(-1,z[c])])
            rhs = -p*s
        else:
            terms.append((-p, z[c]))
            rhs = -p*p*s
        print(f" cusp{c}:", linear(terms), "<=", rhs)
    print("Bounds\n e1 = 1")
    extra = dict(fixed)
    assert len(extra) == len(fixed) and 1 not in extra
    for a in range(2, h + 1):
        if a in extra:
            print(f" {e[a]} = {extra[a]}")
        else:
            print(f" {0 if s >= 3 else -1} <= {e[a]} <= 1")
    for c in indices:
        print(f" 0 <= {z[c]} <= {positive_mass}")
    print(f" -{h} <= multiple <= {h}")
    print("Generals\n", " ".join(e.values()), "multiple\nEnd")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("s", type=int, choices=range(5))
    parser.add_argument("--gap", type=int)
    parser.add_argument("--endpoint", type=int, choices=(-1,1))
    parser.add_argument("--reduced", action="store_true")
    args = parser.parse_args()
    write_model(args.s, args.gap, args.endpoint, args.reduced)
