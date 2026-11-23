# A counterexample to Conjecture 3.5 of Derickx and Stoll

Let $\mu_{\mathbf Q}(p)$ be the least degree of a nonconstant rational modular unit on $X_1(p)$. At $p=191$ one has
$$
\operatorname{gon}_{\mathbf Q} X_1(191)\ \leq\ 380\ <\ \mu_{\mathbf Q}(191)\ \leq\ 478 .
$$
This disproves Conjecture 3.5 of Derickx and Stoll. The strict inequality is computer-assisted. [Conjecture 3.5](https://arxiv.org/html/2505.14109v2#S3)

The map $X_1(191)\to X_0(191)$ has degree $95$. The Fricke quotient $X_0^+(191)$ has genus two, so its canonical map has degree two over $\mathbf Q$. The composite of the three maps has degree $95\cdot2\cdot2=380$. For the other upper bound, the modular unit $E_4E_7/(E_1E_8)$ has degree $478$ by the cusp formula below. [Najman and Orlić, Theorem 1.2](https://arxiv.org/pdf/2207.11650#page=2) · [Elkies, explicit hyperelliptic map](https://people.math.harvard.edu/~elkies/xisog.pdf) · [Derickx and van Hoeij’s table](https://www.math.fsu.edu/~hoeij/files/X1N/gonality)

## Cuspidal orders

Put $p=191$, $h=(p-1)/2=95$ and $D=380$. Index the diamond group by $1,\ldots,95$, identifying $a$ and $-a$ modulo $p$. By Streng's Theorem 1.2, every rational modular unit is, up to a constant, a product of Siegel functions with integer exponents $e_a$ satisfying
$$
\sum_a e_a=12s,\qquad \sum_a a^2 e_a=0\pmod p .
$$
There are $95$ rational cusps and one closed cusp of degree $95$. At the rational cusps the vector of orders is $v=Ae$, where
$$
A_{ca}=\frac p2 B_2\!\left(\left\{\frac{ca}{p}\right\}\right),
\qquad B_2(x)=x^2-x+\frac16 .
$$
This follows from the leading exponent in Streng's equation (1.4) at the cusp with parameter $q^{1/p}$, followed by diamond translation. Each column of $A$ sums to $-h/12$, so the order at the remaining cusp is $s$. Replacing the unit by its inverse, we can assume $s\geq0$. Its degree is $95s+\sum_c\max(v_c,0)$, so $s\leq4$ for any unit of degree at most $380$. [Streng, Theorem 1.2](https://arxiv.org/pdf/1503.08127#page=2)

Let $M=A^{-1}$. For a unit of degree at most $D$, the positive and negative parts of $v$ have mass at most $D-hs$ and $D$, respectively. Their masses differ by $hs$. Adding equal mass at the same coordinate to both parts gives, for any row vector $r$,
$$
r v\ \leq\ (D-hs)\max_i r_i-D\min_i r_i .
$$
Applying this inequality to a row of $M$ and its negative gives the following inclusive integer bounds. Every row is a permutation of the first.

| $s$ | 0 | 1 | 2 | 3 | 4 |
|---|---:|---:|---:|---:|---:|
| Lower bound for $e_a$ | $-2$ | $-1$ | $-1$ | $0$ | $0$ |
| Upper bound for $e_a$ | $2$ | $1$ | $1$ | $1$ | $1$ |

For $s=0$ the same inequality gives $16e_1-e_2<33$ and $4e_1+e_2<8$. If $e_1=2$, these imply $-1<e_2<0$. Diamond translation and inversion therefore exclude every exponent of absolute value two. Thus every exponent belongs to $\{-1,0,1\}$. The script `cusp_bounds.py` verifies these matrix inequalities over $\mathbf Q$. It also gives a separate exclusion for $s=4$.

## Exact finite computation

Write $r_{ca}$ for the integer in $[1,p-1]$ congruent to $ca$ modulo $p$. The rational cusp orders are
$$
v_c=ps+\frac1p\sum_a \frac{r_{ca}(r_{ca}-p)}2\, e_a .
$$
The bound on the degree is equivalent to the existence of real numbers $z_c$ satisfying $z_c\geq0$, $z_c\geq v_c$ and $\sum_c z_c\leq380-95s$. Clearing denominators gives the integer-coefficient model in `unit_model.py`. The congruence is imposed by one additional integer variable $T$, using $\sum_a (a^2\bmod p)e_a=pT$. Since $|e_a|\leq1$, the range $-95\leq T\leq95$ suffices. The program `audit_mps.py` checks every coefficient in the solver input against the Bernoulli formula, together with the variable bounds, integrality markers and symmetry constraints.

Multiplication by $2$ is a cycle of length $95$ on the diamond group. Take gaps in this cyclic order.

- For $s=3,4$, diamond translation puts a positive exponent at $e_1=1$. Each order gives one feasibility problem.
- For $s=1,2$, there are at least $12s$ positive exponents. Their smallest cyclic gap $j$ is at most $\lfloor95/(12s)\rfloor$. A diamond puts two consecutive positives at positions $0,j$. For each distance $t<j$, two positions at distance $t$ cannot both be positive. This gives seven problems for $s=1$ and three for $s=2$.
- For $s=0$, take the smallest gap between **nonzero** exponents. If it is at most seven, a diamond and inversion put $+1$ at position zero and either sign at the other endpoint. For each distance smaller than the gap impose $|e_a|+|e_b|\leq1$, written as four linear inequalities. This gives fourteen problems.
- If $s=0$ and all gaps are at least eight, there are at most eleven nonzero exponents. Their number is even. The program `wide_gaps.cpp` enumerates the cyclic gaps for each possible even support size and all balanced signs up to inversion. It retains the products satisfying the square-sum congruence. The minimum degree it finds is $501$.

The exact solver finds no feasible point in any of the 26 integer-programming cases. For the gap-one balanced case with second sign $-1$, the recorded verification splits by the first nonzero exponent in the list $(3,95,94,61,41,5,8,64)$, and in its all-zero branch uses the second list $(4,63,38,48,32,62,76)$. The all-zero branch of the second list is included. This exhausts the nonconstant rational modular units of degree at most $380$ and proves the stated strict inequality.

## Reproduction

The calculation uses SageMath for rational matrix arithmetic, a C++17 compiler, GLPK to write fixed MPS, and SCIP 10.1.0 with SoPlex 8.1.0 in exact mode. Exact mode is enabled before the input is read. Each solver run in `logs` returned infeasible in exact mode. The case list and node counts are in `checks.tsv`. The program `check_models.py` reports a time limit as unresolved and returns a failure status if any case is unresolved. [SCIP exact-mode documentation](https://scipopt.org/doc/html/EXACT.php)

```sh
sage -python cusp_bounds.py
c++ -O3 -std=c++17 wide_gaps.cpp -o wide_gaps
./wide_gaps
python check_models.py 4
python check_models.py 3
python check_models.py 2 --gaps 1 2 3 --workers 3 --seconds 600
python check_models.py 1 --gaps 1 2 3 4 5 6 7 --workers 4 --seconds 900
python check_models.py 0 --gaps 2 3 4 5 6 7 --nonzero-gap --workers 4 --seconds 900
python check_models.py 0 --gaps 1 --nonzero-gap --endpoint 1 --seconds 900
```

The following command reproduces the last case with both lists in a single partition. The cases with a nonzero entry in the first list impose no restriction on entries in the second list. A larger wall-clock limit may be needed on slower hardware.

```sh
python check_models.py 0 --gaps 1 --nonzero-gap --endpoint -1 --prefix 3 95 94 61 41 5 8 64 4 63 38 48 32 62 76 --workers 4 --seconds 1200
```

As a separate check, `balanced_units.cpp` finds the minimum degree of a nonconstant balanced product with at most five factors on each side, allowing repetitions. For five factors it checks $210388301832$ normalized pairs and obtains $478$. The proof above only needs the short enumeration in `wide_gaps.cpp`.
