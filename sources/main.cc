#include "qc.h"

int main()
{

    int n = 3;
    qc::State psi = qc::basis(n, 0); // |000>
    
    qc::C H[2][2]; qc::gate_H(H);
    qc::apply_1q(H, psi, /*target=*/2); // MSB に H → (|000>+|100>)/√2
    
    qc::pretty_print(psi, n, /*max_terms=*/8, /*cutoff=*/1e-9,
                     /*precision=*/6, /*show_prob=*/true, /*show_phase=*/true);

    qc::C X[2][2]; qc::gate_X(X);
    qc::apply_controlled_1q(X, psi, /*control=*/2, /*target=*/0); // CNOT(2->0)
    
    qc::pretty_print(psi, n, /*max_terms=*/8, /*cutoff=*/1e-9,
                 /*precision=*/6, /*show_prob=*/true, /*show_phase=*/true);

    qc::C U4[4][4]; qc::gate_CNOT(U4);
    qc::apply_2q(U4, psi, /*qA=*/2, /*qB=*/0); // 任意の 2qubit を(2,0)へ

    qc::pretty_print(psi, n, /*max_terms=*/8, /*cutoff=*/1e-9,
                 /*precision=*/6, /*show_prob=*/true, /*show_phase=*/true);
    
}
