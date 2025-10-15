#include "surface_code.h"
#include <numbers>  // if you need std::numbers::pi

namespace qc::surface {

void reset_to_zero(State& psi, int q){
    int m = measure_qubit_Z(psi, q);
    if (m == 1) { C X[2][2]; gate_X(X); apply_1q(X, psi, q); }
}

// Non-destructive: from |0>^9, just apply H to make |+>^9.
void prepare_all_plus_unitary(State& psi) {
    C Hm[2][2]; gate_H(Hm);
    for (int d = 0; d < 9; ++d) apply_1q(Hm, psi, d);
}

// Destructive: Z-measure + X reset → |0> then H → |+>.
void prepare_all_plus_fresh(State& psi) {
    C Hm[2][2]; gate_H(Hm);
    for (int d = 0; d < 9; ++d) {
        reset_to_zero(psi, d);
        apply_1q(Hm, psi, d);
    }
}

// Z round (CNOT data -> anc, then Z-measure on anc)
std::array<int,2> z_round(State& psi) {
    std::array<int,2> syn{0,0};
    C X1[2][2]; gate_X(X1);
    for (int k=0;k<2;++k){
        const int anc = Z_ANC[k];
        reset_to_zero(psi, anc);                     // anc = |0>
        for (int j=0; j<4; ++j) {
            const int d = Z_CHECKS[k][j];
            apply_controlled_1q(X1, psi, /*control=*/d, /*target=*/anc);
        }
        syn[k] = measure_qubit_Z(psi, anc);
    }
    return syn;
}

// X round (anc in |+>, CNOT anc -> data, H, then Z-measure on anc)
std::array<int,2> x_round(State& psi) {
    std::array<int,2> syn{0,0};
    C Hm[2][2]; gate_H(Hm);
    C X1[2][2]; gate_X(X1);
    for (int k=0;k<2;++k){
        const int anc = X_ANC[k];
        reset_to_zero(psi, anc);                     // anc = |0>
        apply_1q(Hm, psi, anc);                      // anc → |+>
        for (int j=0; j<4; ++j) {
            const int d = X_CHECKS[k][j];
            apply_controlled_1q(X1, psi, /*control=*/anc, /*target=*/d);
        }
        apply_1q(Hm, psi, anc);                      // X-measure via H + Z
        syn[k] = measure_qubit_Z(psi, anc);
    }
    return syn;
}

} // namespace qc::surface
