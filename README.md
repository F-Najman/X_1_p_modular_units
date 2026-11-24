# Gonality and modular units on X_1(p)

Code, solver inputs and logs for the paper

> F. Najman, *The Q-gonality of X_1(p) is not always attained by a modular unit*.

The paper proves that gon_Q X_1(191) <= 380, while every non-constant modular unit on
X_1(191) defined over Q has degree at least 381. So the Q-gonality of X_1(p) is not
always attained by a modular unit: this disproves Conjecture 3.5 of Derickx and Stoll
([arXiv:2505.14109v2](https://arxiv.org/abs/2505.14109v2)) for p = 191 and answers
Question 1 of Derickx and van Hoeij (J. Algebra 417 (2014)) negatively for N = 191.
The paper also recalls gon_Q X_1(191) >= 304 (Proposition 3.3, which is Corollary 4.3 of
[hexagon:2611.00884v1](https://hexagonmath.org/2611.00884v1); no computations), so
304 <= gon_Q X_1(191) <= 380 < 381 <= mu_Q(191) <= 478.

All commands below are run from the repository root.

## Requirements

| Tool | Version used | Used for |
|---|---|---|
| SageMath (its `python3`) | 10.7 | exact linear algebra, model generation, HiGHS via scipy 1.17.1 |
| GLPK (`glpsol`) | 5.0 | converting the LP files to fixed MPS |
| SCIP with SoPlex, built with exact solving | 10.1.0 / 8.1.0 | the 56 integer linear programs |
| Magma | V2.30-1 | genera and models of X_0^+(p) |
| C and C++ compilers | gcc 15.2 | the wide-gap enumerations |

The path to the SCIP binary is passed in the environment variable `SCIP`.

## Statement by statement

| Statement in the paper | What is computed | Script | Log |
|---|---|---|---|
| Proposition 3.1 | X_0^+(191) has genus 2 (Claim 5), model y^2+(x^3+x+1)y = x^2-2x (Claim 6) | `code/X0plus.m` | `logs/checks/X0plus.txt` |
| Remark 3.2 | X_0^+(167) has genus 2 (Claim 5), X_0^+(239) has genus 3 and is not hyperelliptic (Claims 5, 7) | `code/X0plus.m` | `logs/checks/X0plus.txt` |
| Introduction | deg H_4H_7/(H_1H_8) = 478 (Claim 1) | `code/independent_checks.py` | `logs/checks/independent_checks.txt` |
| Proposition 5.2, Table 1 | A is invertible, the bounds on the exponents, exclusion of e_k = 2 for s = 0 | `code/cusp_bounds.py`, checked by `code/independent_checks.py` (Claims 1-3) | `logs/rerun_2026-11-23/cusp_bounds.txt`, `logs/checks/independent_checks.txt` |
| Section 6, first paragraph | 2 generates (Z/191Z)^*/{+-1} (Claim 4) | `code/independent_checks.py` | `logs/checks/independent_checks.txt` |
| Proposition 6.2 | the wide-gap case: 3 312 502 exponent vectors, smallest degree 501 | `code/wide_gaps_check.c`, `code/wide_gaps.cpp` | `logs/checks/wide_gaps_check.txt`, `logs/rerun_2026-11-23/wide_gaps.txt` |
| Proposition 6.4, Table 2 | the 56 integer linear programs are infeasible (exact SCIP) | `code/unit_model.py` writes them, `models/*.mps`, `code/rerun_exact.py` solves them | `logs/rerun_2026-11-23/*.log`, `logs/rerun_2026-11-23/summary.tsv` |

## Further checks

These checks are not needed for the proof. Section 8 of the paper mentions the HiGHS rebuild,
the two enumerations and the known positive case.

| Check | Script | Log |
|---|---|---|
| regeneration and row-by-row audit of the 56 inputs | `code/regenerate_models.py`, `code/audit_mps.py` | `logs/checks/regenerate_models.txt` |
| original run of the 56 programs | `code/check_models.py` | `logs/original_run/` |
| the 56 programs rebuilt independently and solved with HiGHS (floating point) | `code/highs_crosscheck.py` | `logs/checks/highs.tsv` |
| independent recomputation of Proposition 5.2 | `code/independent_checks.py` | `logs/checks/independent_checks.txt` |
| two independent enumerations for Proposition 6.2 (counts agree up to the factor 10, 8, 6, 4 of cyclic shifts, except for the two supports with 10 non-zero exponents and gaps alternately 8, 11 and 9, 10, which are invariant under the shift by 19) | `code/wide_gaps_check.c`, `code/wide_gaps.cpp` | as above |
| known positive: with 478 in place of 380, exact SCIP finds (H_4H_7/(H_1H_8))_4 | `code/known_positive.py` | `logs/checks/known_positive.txt` |
| for primes 11 <= p <= 83, the smallest degree of a non-constant unit with all \|e_k\| <= 3 equals the value of Derickx and van Hoeij (gon_Q X_1(p) for p <= 37, their best modular unit for 41 <= p <= 83). This tests the model; for p >= 41 it does not determine mu_Q(p) | `code/validate_small_primes.py` | `logs/checks/validate_small_primes.txt` |

`checks.tsv` lists the 56 programs of Proposition 6.4 (the cusp order s, the gap j, the sign
sigma, the exponents fixed in the split of P(0,1,-1)), the number of nodes of the original run
and the SHA-256 of each stored input in `models/`.

## Reproducing the computations

```sh
python3 code/cusp_bounds.py                     # Proposition 5.2
python3 code/independent_checks.py              # Claims 1-4
cc -O2 -o wide_gaps_check code/wide_gaps_check.c && ./wide_gaps_check            # Proposition 6.2
c++ -O3 -std=c++17 code/wide_gaps.cpp -o wide_gaps && ./wide_gaps                 # Proposition 6.2, second program
python3 code/regenerate_models.py               # regenerate and audit the 56 inputs
SCIP=/path/to/scip python3 code/rerun_exact.py logs/rerun --workers 8            # Proposition 6.4 (about 48 CPU minutes)
SCIP=/path/to/scip python3 code/known_positive.py
python3 code/highs_crosscheck.py logs/checks/highs.tsv 8
python3 code/validate_small_primes.py 83 3
magma -b code/X0plus.m                          # Proposition 3.1, Remark 3.2
```

`code/rerun_exact.py` enables exact solving before reading each problem and counts a program as
excluded only if the log contains `exact/enable = TRUE` and the status
`problem is solved [infeasible]`. The original driver `code/check_models.py` writes the models
with `code/unit_model.py` and solves them in the same way.

## Notes

- Nine of the stored inputs (`order1_gap{1,3,4,5,6,7}`, `order2_gap{1,2,3}`) were written by an
  earlier version of `code/unit_model.py` which left the endpoint row unnamed, so `glpsol` called
  it `R0000004`. They are identical to the regenerated inputs up to this row name, and both
  versions pass the audit (`logs/checks/regenerate_models.txt`).
- `code/balanced_units.cpp` is a supplementary search from the original run (balanced products of
  at most five Siegel functions on each side). It is not used in the paper and was not rerun.
- The computations were run on the Mordell server, Department of Mathematics, University of
  Zagreb (AMD EPYC 9175F, 16 cores, 384 GB RAM).
