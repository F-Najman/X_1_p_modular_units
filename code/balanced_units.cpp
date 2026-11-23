// Balanced products of four or five Siegel functions on X_1(191), with repetition.
// Compile with c++ -O3 -march=native -fopenmp -std=c++17 balanced_units.cpp -o balanced_units
// Run with an argument of 4 or 5. OpenMP and AVX2 are optional.
// Equal-length products with fewer factors are included by padding both sides.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#ifdef __AVX2__
#include <immintrin.h>
#endif

constexpr int p = 191;
constexpr int h = (p - 1) / 2;
constexpr int width = 96;
using Orders = std::array<unsigned char, width>;

int twice_p_order(int a, int c) {
    int r = a * c % p;
    return r * (r - p);  // The omitted constant cancels between products.
}

int total_variation(const Orders& a, const Orders& b) {
#ifdef __AVX2__
    auto sum = _mm256_setzero_si256();
    for (int j = 0; j < width; j += 32) {
        auto x = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(a.data() + j));
        auto y = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(b.data() + j));
        sum = _mm256_add_epi64(sum, _mm256_sad_epu8(x, y));
    }
    alignas(32) std::uint64_t lanes[4];
    _mm256_store_si256(reinterpret_cast<__m256i*>(lanes), sum);
    return static_cast<int>(lanes[0] + lanes[1] + lanes[2] + lanes[3]);
#else
    int sum = 0;
    for (int j = 0; j < width; ++j) sum += std::abs(int(a[j]) - int(b[j]));
    return sum;
#endif
}

template<int factors> int search() {
    using Product = std::array<int, factors>;
    // Since -1 is a nonsquare, simultaneous multiplication of all indices
    // by a diamond takes every common square sum to 0, 1, or -1 modulo p.
    std::array<std::vector<Product>, 3> groups;
    Product indices{};
    auto enumerate = [&](auto& self, int depth, int first, int square_sum) -> void {
        if (depth == factors) {
            int residue = square_sum % p;
            if (residue == 0 || residue == 1 || residue == p-1)
                groups[residue == p-1 ? 2 : residue].push_back(indices);
            return;
        }
        for (int a = first; a <= h; ++a) {
            indices[depth] = a;
            self(self, depth + 1, a, square_sum + a*a);
        }
    };
    enumerate(enumerate, 0, 1, 0);
    int minimum = std::numeric_limits<int>::max();
    Product numerator{}, denominator{};
    std::uint64_t pairs = 0;
    for (const auto& group : groups) {
        std::vector<Orders> signatures;
        for (const auto& product : group) {
            Orders orders{};
            int sum = 0;
            for (int c = 1; c <= h; ++c) {
                int delta = 0;
                for (int i = 0; i < factors; ++i)
                    delta += twice_p_order(product[i], c) - twice_p_order(group[0][i], c);
                assert(delta % (2*p) == 0);
                int order = delta / (2*p);
                assert(-128 <= order && order < 128);
                orders[c-1] = static_cast<unsigned char>(order + 128);
                sum += order;
            }
            assert(sum == 0);
            signatures.push_back(orders);
        }
        #pragma omp parallel
        {
            int best = std::numeric_limits<int>::max();
            size_t best_i = 0, best_j = 0;
            #pragma omp for schedule(dynamic, 64)
            for (size_t i = 0; i < group.size(); ++i) {
                for (size_t j = 0; j < i; ++j) {
                    int norm = total_variation(signatures[i], signatures[j]);
                    if (norm < best) {
                        best = norm;
                        best_i = i;
                        best_j = j;
                    }
                }
            }
            #pragma omp critical
            if (best < minimum) {
                minimum = best;
                numerator = group[best_i];
                denominator = group[best_j];
            }
        }
        pairs += std::uint64_t(group.size()) * (group.size()-1) / 2;
        std::cout << "Normalized bucket size " << group.size() << std::endl;
    }
    assert(minimum > 0 && minimum % 2 == 0);
    std::cout << "Pairs checked " << pairs << "\nMinimum degree " << minimum/2;
    std::cout << "\nNumerator indices";
    for (int a : numerator) std::cout << " " << a;
    std::cout << "\nDenominator indices";
    for (int a : denominator) std::cout << " " << a;
    std::cout << '\n';
    return 0;
}

int main(int argc, char** argv) {
    if (argc == 2 && std::string(argv[1]) == "4") return search<4>();
    if (argc == 2 && std::string(argv[1]) == "5") return search<5>();
    std::cerr << "Expected a factor count of 4 or 5\n";
    return 1;
}
