// Genus of X_0^+(p) = X_0(p)/w_p for p = 167, 191, 239, and a model of X_0^+(191).
// Run from the repository root:  magma -b code/X0plus.m
SetColumns(0);

// Claim 5: X_0^+(191) has genus 2 (and X_0^+(167), X_0^+(239) have genus 2 and 3).
for p in [167, 191, 239] do
    S := CuspForms(p, 2);
    W := AtkinLehnerOperator(S, p);
    g := Dimension(Kernel(W - 1));
    printf "genus of X_0^+(%o) = %o, genus of X_0(%o) = %o\n", p, g, p, Dimension(S);
    assert g eq (p eq 239 select 3 else 2);
end for;
print "Claim 5 successfully verified";

// Claim 6: a model of X_0^+(191), from the q-expansions of a basis of the w_191-invariant
// cusp forms: x = f1/f2, y = q (dx/dq) / f2 satisfy y^2 = F(x) with F of degree 6.
p := 191;
S := CuspForms(p, 2);
W := AtkinLehnerOperator(S, p);
prec := 200;
V := Kernel(W - 1);
B := [&+[ (Basis(V)[i])[j] * qExpansion(Basis(S)[j], prec) : j in [1..Dimension(S)]] : i in [1..2]];
R<q> := Parent(B[1]);
f1 := B[1]; f2 := B[2];
x := f1 / f2;
y := q * Derivative(x) / f2;
P<t> := PolynomialRing(Rationals());
M := Matrix(Rationals(), [ [Coefficient(x^i, n) : n in [-12..prec-20]] : i in [0..6] ] cat
                         [ [Coefficient(y^2, n) : n in [-12..prec-20]] ]);
K := Kernel(M);
assert Dimension(K) eq 1;
k := Basis(K)[1];
F := -&+[ k[i+1] * t^i : i in [0..6] ] / k[8];
C := HyperellipticCurve(F);
C := ReducedMinimalWeierstrassModel(C);
printf "X_0^+(191): %o\n", C;
assert Genus(C) eq 2;
print "Claim 6 successfully verified";

// Claim 7: X_0^+(239) is not hyperelliptic: the 6 products of a basis f1, f2, f3 of the
// w_239-invariant cusp forms are linearly independent. A quadratic relation would be a
// weight 4 form on Gamma_0(239), so it suffices to check this modulo q^81 (Sturm bound
// 4 * 240 / 12 = 80). Hence the canonical model of X_0^+(239) is a plane quartic over Q,
// and projecting from the rational cusp gives a map of degree 3 to P^1 over Q.
p := 239;
S := CuspForms(p, 2);
W := AtkinLehnerOperator(S, p);
V := Kernel(W - 1);
prec := 100;
B := [&+[ (Basis(V)[i])[j] * qExpansion(Basis(S)[j], prec) : j in [1..Dimension(S)]] : i in [1..3]];
prods := [B[i]*B[j] : i in [j..3], j in [1..3]];
M := Matrix(Rationals(), [[Coefficient(f, n) : n in [1..80]] : f in prods]);
assert Rank(M) eq 6;
print "Claim 7 successfully verified";
quit;
