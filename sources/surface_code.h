#pragma once
#include "qc.h"
#include <array>
#include <vector>
#include <cassert>

namespace qc::surface {

// 2d index -> 1d index
inline int data_idx(int i, int j, int d) { return i * d + j; }

// Two Z-faces (TL, BR) and two X-faces (TR, BL) — checkerboard, all weight-4.
struct SurfaceCode {
    int d; // code-length
    int n_data; // = d*d
    std::vector<int> z_anc;
    std::vector<int> x_anc;

    std::vector<std::array<int,4>> z_checks; // z_checks[k] is pair with z_anc[k]
    std::vector<std::array<int,4>> x_checks; // x_checks[k] is pair with x_anc[k]

    int n_qubits() const { return n_data + (int)z_anc.size() + (int)x_anc.size(); }

};

SurfaceCode build_surface_code(int d);

// Measure in Z and (if needed) apply X to force |0>.
void reset_to_zero(State& psi, int q);

// Prepare |+>^9 non-destructively (assumes |0>^9 → just H on data 0..8).
void prepare_all_plus_unitary(State& psi, const SurfaceCode& sc);

// Prepare |+>^9 destructively (Z-measure + X reset on each data, then H).
// WARNING: This erases any pre-existing errors/phases on data qubits.
// Use only at the start of an independent run, before injecting errors.
void prepare_all_plus_fresh(State& psi, const SurfaceCode& sc);

// One Z stabilizer round: anc in |0>, CNOT(data -> anc), Z-measure.
std::vector<int> z_round(State& psi, const SurfaceCode& sc);

// One X stabilizer round (standard): anc in |+>, CNOT(anc -> data), H, Z-measure.
std::vector<int> x_round(State& psi, const SurfaceCode& sc);

} // namespace qc::surface
