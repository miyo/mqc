#include "qc.h"
#include "surface_code.h"

#include <array>
#include <vector>

namespace qc::surface {

void reset_to_zero(State& psi, int q){
    int m = measure_qubit_Z(psi, q);
    if (m == 1) { C X[2][2]; gate_X(X); apply_1q(X, psi, q); }
}

// Non-destructive: from |0>^9, just apply H to make |+>^9.
void prepare_all_plus_unitary(State& psi, const SurfaceCode& sc) {
    C Hm[2][2]; gate_H(Hm);
    for (int d = 0; d < sc.n_data; ++d) apply_1q(Hm, psi, d);
}

// Destructive: Z-measure + X reset → |0> then H → |+>.
void prepare_all_plus_fresh(State& psi, const SurfaceCode& sc) {
    C Hm[2][2]; gate_H(Hm);
    for (int d = 0; d < sc.n_data; ++d) {
        reset_to_zero(psi, d);
        apply_1q(Hm, psi, d);
    }
}

// Z round (CNOT data -> anc, then Z-measure on anc)
std::vector<int> z_round(State& psi, const SurfaceCode& sc) {
    std::vector<int> syn(sc.z_anc.size(), 0);
    C X1[2][2]; gate_X(X1);
    for (size_t k = 0; k < sc.z_anc.size(); ++k){
        const int anc = sc.z_anc[k];
        reset_to_zero(psi, anc); // anc = |0>
        const auto& nb = sc.z_checks[k]; // target data qubits
        for (int t = 0; t < 4; ++t) {
            const int dqb = nb[t];
            apply_controlled_1q(X1, psi, /*control=*/dqb, /*target=*/anc);
        }
        syn[k] = measure_qubit_Z(psi, anc);
    }
    return syn;
}

// X round (anc in |+>, CNOT anc -> data, H, then Z-measure on anc)
std::vector<int> x_round(State& psi, const SurfaceCode& sc) {
    std::vector<int> syn(sc.x_anc.size(), 0);
    C Hm[2][2]; gate_H(Hm);
    C X1[2][2]; gate_X(X1);
    for (size_t k = 0; k < sc.x_anc.size(); ++k){
        const int anc = sc.x_anc[k];
        reset_to_zero(psi, anc); // anc = |0>
        apply_1q(Hm, psi, anc); // anc → |+>
        const auto& nb = sc.x_checks[k]; // target data qubits
        for (int t = 0; t < 4; ++t) {
            const int dqb = nb[t];
            apply_controlled_1q(X1, psi, /*control=*/anc, /*target=*/dqb);
        }
        apply_1q(Hm, psi, anc); // X-measure via H + Z
        syn[k] = measure_qubit_Z(psi, anc);
    }
    return syn;
}

SurfaceCode build_surface_code(int d) {
    assert(d >= 3 && (d % 2 == 1));
    SurfaceCode sc;
    sc.d = d;
    sc.n_data = d * d;

    std::vector<std::array<int,4>> z_checks, x_checks;

    // (d-1)×(d-1) のプラケット格子をチェッカーボード色分け
    for (int i = 0; i < d - 1; ++i) {
        for (int j = 0; j < d - 1; ++j) {
            std::array<int,4> nbrs{
                data_idx(i,   j,   d),
                data_idx(i+1, j,   d),
                data_idx(i,   j+1, d),
                data_idx(i+1, j+1, d)
            };
            if ( ((i + j) & 1) == 0 ) {
                z_checks.push_back(nbrs);
            } else {
                x_checks.push_back(nbrs);
            }
        }
    }

    // 物理インデックスの付与：data [0..d*d-1], 次に Z anc, 次に X anc
    int next = sc.n_data;
    sc.z_anc.resize(z_checks.size());
    for (size_t k = 0; k < z_checks.size(); ++k) sc.z_anc[k] = next++;
    sc.x_anc.resize(x_checks.size());
    for (size_t k = 0; k < x_checks.size(); ++k) sc.x_anc[k] = next++;

    sc.z_checks = std::move(z_checks);
    sc.x_checks = std::move(x_checks);
    return sc;
}

} // namespace qc::surface
