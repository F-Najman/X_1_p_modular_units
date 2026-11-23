"""Check the exact solver input, including every row, bound and integer marker."""
from collections import defaultdict
from fractions import Fraction


def audit_mps(path, s, gap, endpoint, reduced=False, fixed=()):
    row_types, rhs, bounds = {}, {}, {}
    by_row = defaultdict(dict)
    integral, columns = set(), set()
    inside_integer = False
    section = None
    for raw in path.read_text().splitlines():
        if not raw or raw[0] == "*":
            continue
        if raw[0] != " ":
            section = raw.split()[0]
            continue
        t = raw.split()
        if section == "ROWS":
            kind, row = t
            assert row not in row_types
            row_types[row] = kind
        elif section == "COLUMNS":
            if t[1] == "'MARKER'":
                assert t[2] in ("'INTORG'", "'INTEND'")
                inside_integer = t[2] == "'INTORG'"
                continue
            var = t[0]
            columns.add(var)
            if inside_integer:
                integral.add(var)
            for name, val in zip(t[1::2], t[2::2]):
                assert name in row_types and var not in by_row[name]
                by_row[name][var] = Fraction(val)
        elif section == "RHS":
            for name, val in zip(t[1::2], t[2::2]):
                assert name in row_types and name not in rhs
                rhs[name] = Fraction(val)
        elif section == "BOUNDS":
            kind, _, var, value = t
            value = Fraction(value)
            limits = bounds.setdefault(var, [0, None])
            assert kind in ("LO", "UP", "FX")
            if kind in ("LO", "FX"):
                limits[0] = value
            if kind in ("UP", "FX"):
                limits[1] = value
        else:
            raise ValueError(f"Unexpected MPS section {section}")

    p, h = 191, 95
    idx = range(1, h+1)
    e, z = {a:f"e{a}" for a in idx}, {c:f"z{c}" for c in idx}
    assert integral == set(e.values()) | {"multiple"}
    assert columns == integral | set(z.values())
    expected_bounds = {e[a]:[0 if s >= 3 else -1, 1] for a in idx}
    expected_bounds[e[1]] = [1, 1]
    for a, v in fixed:
        expected_bounds[e[a]] = [v, v]
    expected_bounds.update({z[c]:[0, 380-95*s] for c in idx})
    expected_bounds["multiple"] = [-95,95]
    assert expected_bounds == bounds

    actual = iter(row_types)

    def expect(kind, coefficients, right=0):
        row = next(actual)
        assert row_types[row] == kind and rhs.get(row, 0) == right
        assert by_row[row] == {v:n for v,n in coefficients.items() if n}

    expect("N", dict.fromkeys(z.values(),1))
    expect("L", dict.fromkeys(z.values(),1), 380-95*s)
    expect("E", dict.fromkeys(e.values(),1), 12*s)
    center = lambda a: a % p if a % p <= 95 else a % p - p
    residues = {a: center(a*a*96) if reduced else a*a % p for a in idx}
    square_row = {e[a]:residues[a] for a in idx}
    square_row["multiple"] = -p
    expect("E", square_row)

    if gap is not None:
        cycle = []
        a = 1
        for _ in idx:
            cycle.append(a)
            a = abs(center(2*a))
        assert len(set(cycle)) == h and a == 1
        expect("E", {e[cycle[gap]]:1}, 1 if endpoint is None else endpoint)
        for shift in range(1, gap):
            for i in range(h):
                left, right = cycle[i], cycle[(i+shift)%h]
                for x,y in ([(1,1)] if endpoint is None else [(1,1),(1,-1),(-1,1),(-1,-1)]):
                    expect("L",{e[left]:x, e[right]:y},1)

    for c in idx:
        row = {}
        factor = center(c*c)
        for a in idx:
            r = a*c % p
            # Start from p B_2(r/p)/2, then subtract its constant p/12.
            bernoulli = (Fraction(r,p)**2 - Fraction(r,p) + Fraction(1,6))*p/2
            coefficient = bernoulli - Fraction(p,12)
            if reduced:
                row[e[a]] = coefficient - Fraction(factor*residues[a],p)
            else:
                row[e[a]] = p*coefficient
            assert row[e[a]].denominator == 1
        row[z[c]] = -1 if reduced else -p
        if reduced:
            row["multiple"] = factor
        expect("L", row, -p*s*(1 if reduced else p))
    assert next(actual, None) is None
