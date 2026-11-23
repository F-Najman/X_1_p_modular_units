/* Independent check of the wide-gap case (s = 0, all exponents in {-1,0,1}, any two
 * nonzero exponents at cyclic distance at least 8 in the order 1, 2, 2^2, ... of the
 * diamond group modulo signs).
 *
 * Written independently of wide_gaps.cpp, with a different enumeration: we do not
 * reduce the support modulo rotations, we only use a diamond to put a nonzero exponent
 * at the position of a = 1, and inversion to make it +1.  Every composition of 95 into
 * k parts of size at least 8 is then enumerated, for k = 2, 4, ..., 10.
 *
 * For every balanced sign vector satisfying sum e_a a^2 = 0 mod 191 we compute the orders
 * v_c at the 95 rational cusps exactly (as integers, after multiplying by 191) and the
 * degree sum_c max(v_c, 0).  We print the minimal degree found.
 *
 * Compile: cc -O2 -o wide_gaps_check wide_gaps_check.c
 */
#include <stdio.h>
#include <stdlib.h>

#define P 191
#define H 95
#define GAP 8
#define BOUND 380

static int idx[H];          /* idx[i] = 2^i mod P folded to [1, H] */
static long pord[H][H];     /* pord[i][c-1] = 2P * (order of H_{idx[i]} at cusp c) - P^2/6 */
static int sq[H];           /* idx[i]^2 mod P */
static int pos[12];
static long minimum = 1L << 40;
static unsigned long long supports = 0, candidates = 0;

static int fold(int x) { x %= P; if (x < 0) x += P; return x <= H ? x : P - x; }

static void check_signs(int k) {
    /* sign of pos[0] is +1; choose (k/2 - 1) further + signs among pos[1..k-1] */
    int m = k - 1;
    for (int mask = 0; mask < (1 << m); mask++) {
        if (__builtin_popcount(mask) != k / 2 - 1) continue;
        int sign[12];
        sign[0] = 1;
        for (int j = 1; j < k; j++) sign[j] = (mask >> (j - 1)) & 1 ? 1 : -1;
        long r = 0;
        for (int j = 0; j < k; j++) r += sign[j] * sq[pos[j]];
        if (r % P) continue;
        candidates++;
        /* v_c * 2P = sum_j sign_j * (r_cj (r_cj - P)) + (P^2/6) * sum_j sign_j, and
         * sum_j sign_j = 0, so v_c = sum_j sign_j * pord[pos_j][c] / (2P). */
        long degree = 0, total = 0;
        for (int c = 0; c < H; c++) {
            long v2p = 0;
            for (int j = 0; j < k; j++) v2p += sign[j] * pord[pos[j]][c];
            if (v2p % (2 * P)) { fprintf(stderr, "non-integral order\n"); exit(1); }
            long v = v2p / (2 * P);
            total += v;
            if (v > 0) degree += v;
        }
        if (total != 0) { fprintf(stderr, "orders do not sum to zero\n"); exit(1); }
        if (degree < minimum) minimum = degree;
    }
}

static void enumerate(int depth, int k) {
    if (depth == k) {
        if (H - pos[k - 1] < GAP) return;  /* cyclic gap back to position 0 */
        supports++;
        check_signs(k);
        return;
    }
    for (int i = pos[depth - 1] + GAP; i <= H - GAP; i++) {
        pos[depth] = i;
        enumerate(depth + 1, k);
    }
}

int main(void) {
    int x = 1;
    for (int i = 0; i < H; i++) {
        idx[i] = fold(x);
        x = 2 * x % P;
    }
    int seen[H + 1] = {0};
    for (int i = 0; i < H; i++) {
        if (seen[idx[i]]) { fprintf(stderr, "2 does not generate (Z/191)^*/{+-1}\n"); return 1; }
        seen[idx[i]] = 1;
    }
    for (int i = 0; i < H; i++) {
        sq[i] = idx[i] * idx[i] % P;
        for (int c = 1; c <= H; c++) {
            int r = idx[i] * c % P;
            pord[i][c - 1] = (long)r * (r - P);
        }
    }
    pos[0] = 0;
    for (int k = 2; GAP * k <= H; k += 2) {
        supports = candidates = 0;
        enumerate(1, k);
        printf("%2d nonzero exponents: %llu supports through position 0, %llu satisfy the congruence\n",
               k, supports, candidates);
    }
    printf("Minimal degree with all nonzero gaps at least %d: %ld\n", GAP, minimum);
    if (minimum <= BOUND) { printf("FAILURE: degree <= %d found\n", BOUND); return 1; }
    printf("WIDE GAP CASE VERIFIED: all degrees > %d\n", BOUND);
    return 0;
}
