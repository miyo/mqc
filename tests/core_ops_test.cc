// tests/core_ops_test.cc
#include "qc.h"
#include <gtest/gtest.h>
#include <array>
#include <vector>
#include <cmath>

using namespace qc;

namespace {
// small helper for comparing complex states
void expect_state_eq(const State& psi,
                     const std::vector<C>& ref,
                     double tol = 1e-12) {
    ASSERT_EQ(psi.size(), ref.size());
    for (size_t i = 0; i < psi.size(); ++i) {
        EXPECT_NEAR(psi[i].real(), ref[i].real(), tol) << "i=" << i << " (real)";
        EXPECT_NEAR(psi[i].imag(), ref[i].imag(), tol) << "i=" << i << " (imag)";
    }
}
} // namespace

// ---------- basis ----------
TEST(CoreOps, Basis_SizeAndOneHot) {
    const int n = 3;
    const std::uint64_t idx = 5; // b0101 (|q2 q1 q0> = |101>)
    State psi = basis(n, idx);

    ASSERT_EQ(psi.size(), (size_t)(1u << n));
    // only idx is 1, others 0
    for (size_t i = 0; i < psi.size(); ++i) {
        if (i == idx) {
            EXPECT_NEAR(std::abs(psi[i]), 1.0, 1e-12);
        } else {
            EXPECT_NEAR(std::abs(psi[i]), 0.0, 1e-12);
        }
    }
}

// ---------- apply_1q ----------
TEST(CoreOps, Apply1Q_XOnTarget) {
    // 3 qubits, start in |000>
    State psi = basis(/*n=*/3, /*index=*/0);

    // X on qubit-1 (LSB = 0)
    C X[2][2]; gate_X(X);
    apply_1q(X, psi, /*target=*/1);

    // Expect |010> (index b010 = 2)
    std::vector<C> ref(8, C{0,0});
    ref[2] = C{1,0};
    expect_state_eq(psi, ref);
}

TEST(CoreOps, Apply1Q_HNormalization) {
    // single qubit |0>
    State psi = basis(1, 0);

    C H[2][2]; gate_H(H);
    apply_1q(H, psi, /*target=*/0);

    const double s = 1.0/std::sqrt(2.0);
    std::vector<C> ref = { C{s,0}, C{s,0} };
    expect_state_eq(psi, ref, 1e-12);
}

// ---------- apply_2q ----------
TEST(CoreOps, Apply2Q_CNOT_MappingOnBasis) {
    // 2 qubits; order is |q1 q0>, LSB = q0
    // Use gate_CNOT which is the standard matrix with control = low bit, target = high bit
    C U4[4][4]; gate_CNOT(U4);

    // Case 1: |00> -> stays |00>
    {
        State psi = basis(2, 0);  // |00>
        apply_2q(U4, psi, /*qA=*/0, /*qB=*/1); // low=0 control, high=1 target
        std::vector<C> ref(4, C{0,0}); ref[0] = C{1,0};
        expect_state_eq(psi, ref);
    }
    // Case 2: |01> -> stays |01> (control=0)
    {
        State psi = basis(2, 1);  // |01>
        apply_2q(U4, psi, 0, 1);
        std::vector<C> ref(4, C{0,0}); ref[1] = C{1,0};
        expect_state_eq(psi, ref);
    }
    // Case 3: |10> -> flips target |11> (control=0)
    {
        State psi = basis(2, 2);  // |10>
        apply_2q(U4, psi, 0, 1);
        std::vector<C> ref(4, C{0,0}); ref[3] = C{1,0};
        expect_state_eq(psi, ref);
    }
    // Case 4: |11> -> flips target => |10>
    {
        State psi = basis(2, 3);  // |11>
        apply_2q(U4, psi, 0, 1);
        std::vector<C> ref(4, C{0,0}); ref[2] = C{1,0};
        expect_state_eq(psi, ref);
    }
}

TEST(CoreOps, Apply2Q_CreateBell) {
    // Bell via: H on control (q0), then CNOT(control=q0 -> target=q1)
    State psi = basis(2, 0); // |00>

    C H[2][2]; gate_H(H);
    apply_1q(H, psi, /*target=*/1); // (|00> + |10>)/sqrt2

    C U4[4][4]; gate_CNOT(U4);       // control = low (q0), target = high (q1)
    apply_2q(U4, psi, 0, 1);         // -> (|00> + |11>)/sqrt2

    const double s = 1.0/std::sqrt(2.0);
    std::vector<C> ref = { C{s,0}, C{0,0}, C{0,0}, C{s,0} };
    expect_state_eq(psi, ref, 1e-12);
}
