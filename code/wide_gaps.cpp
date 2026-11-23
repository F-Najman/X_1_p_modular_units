// Check all balanced ternary exponent vectors whose nonzeros, in the diamond
// order 1,2,2^2,... modulo signs, have cyclic gaps at least eight.
// Compile with c++ -O3 -std=c++17 wide_gaps.cpp -o wide_gaps
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

int main() {
    constexpr int p = 191, h = 95;
    std::array<int,h> index{}, square{};
    std::array<std::array<int,h>,h> order{};
    int power = 1;
    for (int i = 0; i < h; ++i) {
        index[i] = std::min(power, p-power);
        square[i] = index[i]*index[i] % p;
        for (int c = 1; c <= h; ++c) {
            int r = index[i]*c % p;
            order[i][c-1] = r*(r-p)/2;
        }
        power = 2*power % p;
    }
    auto distinct = index;
    std::sort(distinct.begin(), distinct.end());
    for (int i = 0; i < h; ++i) assert(distinct[i] == i+1);
    int minimum = std::numeric_limits<int>::max();
    for (int size = 2; 8*size <= h; size += 2) {
        std::uint64_t supports = 0, modular_products = 0;
        std::vector<int> positions(size,0), gaps(size);
        std::vector<std::vector<int>> signs;
        for (int mask = 1; mask < (1<<size); mask += 2) {
            std::vector<int> sign(size, -1);
            for (int j = 0; j < size; ++j)
                if (mask & (1<<j)) sign[j] = 1;
            if (std::count(sign.begin(), sign.end(), 1) == size/2)
                signs.push_back(sign);
        }
        auto enumerate = [&](auto& self, int depth) -> void {
            if (depth < size) {
                int last = h - 8*(size-depth);
                for (int i = positions[depth-1]+8; i <= last; ++i) {
                    positions[depth] = i;
                    self(self, depth+1);
                }
                return;
            }
            for (int j = 1; j < size; ++j) gaps[j-1] = positions[j]-positions[j-1];
            gaps.back() = h-positions.back();
            // Keep lexicographically smallest rotations of the gap word.
            for (int t = 1; t < size; ++t) {
                int j = 0;
                while (j < size && gaps[(t+j)%size] == gaps[j]) ++j;
                if (j < size && gaps[(t+j)%size] < gaps[j]) return;
            }
            ++supports;
            for (const auto& sign : signs) {
                int residue = 0;
                for (int j = 0; j < size; ++j) residue += sign[j]*square[positions[j]];
                if (residue % p) continue;
                ++modular_products;
                int degree = 0, sum = 0;
                for (int c = 0; c < h; ++c) {
                    int v = 0;
                    for (int j = 0; j < size; ++j) v += sign[j]*order[positions[j]][c];
                    assert(v % p == 0);
                    sum += v/p;
                    degree += std::max(0, v/p);
                }
                assert(sum == 0 && degree > 380);
                minimum = std::min(minimum,degree);
            }
        };
        enumerate(enumerate,1);
        std::cout << size << " nonzeros: " << supports << " supports, "
                  << modular_products << " modular products\n";
    }
    std::cout << "Minimum degree for gaps at least eight: " << minimum << '\n';
}
