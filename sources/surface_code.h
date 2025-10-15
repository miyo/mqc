#pragma once
#include "qc.h"
#include <array>

namespace qc::surface {

// 3x3 data qubits laid out row-major; center = 4.
inline constexpr int d_idx(int r, int c) { return 3*r + c; }

// Two Z-faces (TL, BR) and two X-faces (TR, BL) — checkerboard, all weight-4.
inline constexpr std::array<int,2> Z_ANC{ 9, 10 };
inline constexpr std::array<int,2> X_ANC{ 11, 12 };

// Use fixed-size arrays (4 neighbors per face) to allow constexpr in header.
inline constexpr std::array<std::array<int,4>,2> Z_CHECKS{{
    { d_idx(0,0), d_idx(0,1), d_idx(1,0), d_idx(1,1) }, // TL
    { d_idx(1,1), d_idx(1,2), d_idx(2,1), d_idx(2,2) }  // BR
}};
inline constexpr std::array<std::array<int,4>,2> X_CHECKS{{
    { d_idx(0,1), d_idx(0,2), d_idx(1,1), d_idx(1,2) }, // TR
    { d_idx(1,0), d_idx(2,0), d_idx(1,1), d_idx(2,1) }  // BL
}};

// --- helpers (declarations) ---

// Measure in Z and (if needed) apply X to force |0>.
void reset_to_zero(State& psi, int q);

// Prepare |+>^9 non-destructively (assumes |0>^9 → just H on data 0..8).
void prepare_all_plus_unitary(State& psi);

// Prepare |+>^9 destructively (Z-measure + X reset on each data, then H).
// WARNING: This erases any pre-existing errors/phases on data qubits.
// Use only at the start of an independent run, before injecting errors.
void prepare_all_plus_fresh(State& psi);

// One Z stabilizer round: anc in |0>, CNOT(data -> anc), Z-measure.
std::array<int,2> z_round(State& psi);

// One X stabilizer round (standard): anc in |+>, CNOT(anc -> data), H, Z-measure.
std::array<int,2> x_round(State& psi);

} // namespace qc::surface
