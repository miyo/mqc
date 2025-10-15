#include "qc.h"

#include <vector>
#include <complex>
#include <random>
#include <cmath>
#include <numbers>
#include <cstdint>
#include <algorithm>

namespace qc {

// Random number generator [0,1)
inline double urand() {
    static std::mt19937_64 eng(std::random_device{}());
    static std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(eng);
}

void renormalize(State& psi) {
    double s2 = 0.0; for (auto& a : psi) s2 += std::norm(a);
    if (s2 <= 0.0) return;                 // already zero vector -> skip
    double s = std::sqrt(s2);
    for (auto& a : psi) a /= s;
}

// construct |basis⟩ 
State basis(int n_qubits, std::uint64_t index) {
    const std::size_t N = 1ull << n_qubits;
    State psi(N, C{0,0});
    if (index < N) psi[index] = C{1,0};
    return psi;
}

// Apply the 1-qubit gate U to the target qubit (O(2^n)).
// Bit numbering: LSB = 0. U is a 2×2 row-major matrix.
void apply_1q(const C U[2][2], State& psi, int target) {
    const std::size_t N   = psi.size();
    const std::size_t step  = 1ull << target;   // 0100...
    const std::size_t block = step << 1;        // 1000...

    for (std::size_t base = 0; base < N; base += block) {
        for (std::size_t off = 0; off < step; ++off) {
            const std::size_t i0 = base + off;       // target bit = 0
            const std::size_t i1 = i0 + step;        // target bit = 1
            const C a = psi[i0], b = psi[i1];
            psi[i0] = U[0][0]*a + U[0][1]*b;
            psi[i1] = U[1][0]*a + U[1][1]*b;
        }
    }
}

// Apply an arbitrary 2-qubit gate U4 (4×4) to qubits (qA, qB) (order-agnostic).
void apply_2q(const C U4[4][4], State& psi, int qA, int qB) {
    const int low  = std::min(qA, qB);
    const int high = std::max(qA, qB);
    const std::size_t sL = 1ull << low;
    const std::size_t sH = 1ull << high;
    const std::size_t N  = psi.size();

    for (std::size_t base = 0; base < N; base += (1ull << (high + 1))) {
        for (std::size_t mid = 0; mid < (1ull << high); mid += (1ull << (low + 1))) {
            for (std::size_t off = 0; off < sL; ++off) {
                const std::size_t i00 = base + mid + off;
                const std::size_t i01 = i00 + sL;
                const std::size_t i10 = i00 + sH;
                const std::size_t i11 = i10 + sL;

                const C v00 = psi[i00], v01 = psi[i01], v10 = psi[i10], v11 = psi[i11];
                C w00 = U4[0][0]*v00 + U4[0][1]*v01 + U4[0][2]*v10 + U4[0][3]*v11;
                C w01 = U4[1][0]*v00 + U4[1][1]*v01 + U4[1][2]*v10 + U4[1][3]*v11;
                C w10 = U4[2][0]*v00 + U4[2][1]*v01 + U4[2][2]*v10 + U4[2][3]*v11;
                C w11 = U4[3][0]*v00 + U4[3][1]*v01 + U4[3][2]*v10 + U4[3][3]*v11;

                psi[i00] = w00; psi[i01] = w01; psi[i10] = w10; psi[i11] = w11;
            }
        }
    }
}

static void make_controlled_U(C U4[4][4], const C U[2][2], bool control_is_high) {
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) U4[i][j] = C{0,0};

    if (control_is_high) {
        // In (high, low) ordering, control = high -> block diag(I_low, U_low)
        U4[0][0] = C{1,0}; U4[1][1] = C{1,0};       // Upper-left I2
        U4[2][2] = U[0][0]; U4[2][3] = U[0][1];     // Lower-right U
        U4[3][2] = U[1][0]; U4[3][3] = U[1][1];
    } else {
        // control = low -> apply U on the high qubit when the low bit = 1
        // Matrix row/col order is [00, 01, 10, 11]
        // Place U on the subspace {01, 11} (indices 1 and 3)
        U4[0][0] = C{1,0}; U4[2][2] = C{1,0};       // When low = 0: identity on {00, 10}
        U4[1][1] = U[0][0]; U4[1][3] = U[0][1];     // When low = 1: apply U on {01, 11}
        U4[3][1] = U[1][0]; U4[3][3] = U[1][1];
    }
}

// This applies a controlled 1-qubit gate via apply_2q (control → target).
void apply_controlled_1q(const C U[2][2], State& psi, int control, int target) {
    C U4[4][4];
    make_controlled_U(U4, U, /*control_is_high=*/ control > target);
    apply_2q(U4, psi, control, target);
}

// Z-measurement
std::uint64_t measure_all(State& psi) {
    // soft normalize
    double s2 = 0.0; for (auto& a : psi) s2 += std::norm(a);
    if (s2 > 0.0) {
        double s = std::sqrt(s2);
        for (auto& a : psi) a /= s;
    }

    double r = urand(), cum = 0.0;
    std::uint64_t idx = psi.empty() ? 0 : (psi.size() - 1);
    for (std::uint64_t i = 0; i < psi.size(); ++i) {
        cum += std::norm(psi[i]);
        if (r < cum) { idx = i; break; }
    }
    for (auto& a : psi) a = C{0,0};
    if (!psi.empty()) psi[idx] = C{1,0};
    return idx;
}

// Measure a single qubit in Z basis and collapse the state.
// Returns 0/1. Collapses in-place and renormalizes the kept subspace.
int measure_qubit_Z(State& psi, int target) {
    const std::size_t N    = psi.size();
    const std::size_t step = 1ull << target;
    const std::size_t block = step << 1;

    // Guard: invalid target (e.g., target >= log2(N))
    if (step == 0 || block == 0 || step >= N) return 0;

    // Compute probabilities for target=0 and target=1
    double n0 = 0.0, n1 = 0.0;
    for (std::size_t base = 0; base < N; base += block) {
        for (std::size_t off = 0; off < step; ++off) {
            n0 += std::norm(psi[base + off]);        // target bit = 0
            n1 += std::norm(psi[base + off + step]); // target bit = 1
        }
    }
    const double denom = n0 + n1;
    if (denom <= 0.0) {
        // Degenerate state: leave |...0> by convention
        return 0;
    }
    double p0 = n0 / denom;

    // Snap near 0/1 to be robust against rounding
    constexpr double eps = 1e-6;
    if (p0 <= eps) p0 = 0.0;
    else if (p0 >= 1.0 - eps) p0 = 1.0;

    // Sample outcome
    const double r = urand();
    const int outcome = (p0 == 0.0) ? 1 :
                        (p0 == 1.0) ? 0 :
                        (r < p0 ? 0 : 1);

    // Collapse and renormalize only the kept half
    const double keep_norm = (outcome == 0) ? n0 : n1;
    const double inv = (keep_norm > 0.0) ? 1.0 / std::sqrt(keep_norm) : 0.0;

    if (outcome == 0) {
        for (std::size_t base = 0; base < N; base += block) {
            // keep i0, scale
            for (std::size_t off = 0; off < step; ++off)
                psi[base + off] *= inv;
            // zero i1
            for (std::size_t off = 0; off < step; ++off)
                psi[base + off + step] = C{0,0};
        }
    } else {
        for (std::size_t base = 0; base < N; base += block) {
            // zero i0
            for (std::size_t off = 0; off < step; ++off)
                psi[base + off] = C{0,0};
            // keep i1, scale
            for (std::size_t off = 0; off < step; ++off)
                psi[base + off + step] *= inv;
        }
    }
    return outcome;
}

}
