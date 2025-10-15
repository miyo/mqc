// tests/measure_and_ctrl_test.cc
#include "qc.h"
#include <gtest/gtest.h>
#include <vector>
#include <cmath>

using namespace qc;

namespace {
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

// ------------------------------------------------------------
// measure_qubit_Z
// ------------------------------------------------------------

TEST(MeasureZ, CollapsesToZeroOnKet0) {
    State psi = basis(/*n=*/1, /*index=*/0);  // |0>
    int m = measure_qubit_Z(psi, /*target=*/0);
    EXPECT_EQ(m, 0);
    std::vector<C> ref = { C{1,0}, C{0,0} };
    expect_state_eq(psi, ref);
}

TEST(MeasureZ, CollapsesToOneOnKet1) {
    State psi = basis(/*n=*/1, /*index=*/1);  // |1>
    int m = measure_qubit_Z(psi, /*target=*/0);
    EXPECT_EQ(m, 1);
    std::vector<C> ref = { C{0,0}, C{1,0} };
    expect_state_eq(psi, ref);
}

// snap（p0 ≈ 1）を強めに検証：p0 = 1 - 1e-12 でも outcome=0 に“決め打ち”
TEST(MeasureZ, SnapNearOne) {
    State psi(2);
    const double delta = 1e-12;
    psi[0] = C{ std::sqrt(1.0 - delta), 0.0 };  // |0> 振幅
    psi[1] = C{ std::sqrt(delta),        0.0 };  // |1> 振幅
    // ここで p0 = 1 - 1e-12。実装側の eps=1e-6 程度なら snap で 0 を返す想定
    int m = measure_qubit_Z(psi, /*target=*/0);
    EXPECT_EQ(m, 0);
    // 完全に |0> へ射影・正規化されているはず
    std::vector<C> ref = { C{1,0}, C{0,0} };
    expect_state_eq(psi, ref, 1e-12);
}

// マルチキュービットの射影挙動：|01> を target=q0 で測ると 1 を返し、
// |01> 以外の振幅が 0 になり、全体が正規化される
TEST(MeasureZ, TwoQubitCollapseOnTargetLSB) {
    // |q1 q0> で |01>
    State psi = basis(/*n=*/2, /*index=*/1);
    int m = measure_qubit_Z(psi, /*target=*/0);  // 測るのは LSB(q0)
    EXPECT_EQ(m, 1);
    std::vector<C> ref(4, C{0,0});
    ref[1] = C{1,0};  // そのまま |01> に残る
    expect_state_eq(psi, ref);
}

// ------------------------------------------------------------
// apply_controlled_1q
// ------------------------------------------------------------

// U=X を与えたときの “CNOT(control -> target)” としての動作を
// control/target の大小（low/high）に関わらず検証

TEST(Controlled1Q, CNOT_ControlHigh_TargetLow) {
    // q1 を control, q0 を target
    // 期待: |00>->|00>, |01>->|01>, |10>->|11>, |11>->|10>
    C Xg[2][2]; gate_X(Xg);

    auto check = [&](std::uint64_t in, std::uint64_t out) {
        State psi = basis(2, in);
        apply_controlled_1q(Xg, psi, /*control=*/1, /*target=*/0);
        State ref = basis(2, out);
        expect_state_eq(psi, ref);
    };
    check(0b00, 0b00);
    check(0b01, 0b01);
    check(0b10, 0b11);
    check(0b11, 0b10);
}

TEST(Controlled1Q, CNOT_ControlLow_TargetHigh) {
    // q0 を control, q1 を target
    // 期待: |00>->|00>, |01>->|11>, |10>->|10>, |11>->|01>
    C Xg[2][2]; gate_X(Xg);

    auto check = [&](std::uint64_t in, std::uint64_t out) {
        State psi = basis(2, in);
        apply_controlled_1q(Xg, psi, /*control=*/0, /*target=*/1);
        State ref = basis(2, out);
        expect_state_eq(psi, ref);
    };
    check(0b00, 0b00);
    check(0b01, 0b11);
    check(0b10, 0b10);
    check(0b11, 0b01);
}

// U=H を渡した場合：control=1 のときのみ target に H がかかる
TEST(Controlled1Q, ControlledH_ActsOnlyWhenControlOne) {
    C Hg[2][2]; gate_H(Hg);

    // 入力 |10>（q1=1, q0=0）。control=q1, target=q0 に H を条件適用。
    // 結果は |1> ⊗ H|0> = |1> ⊗ (|0>+|1>)/√2
    {
        State psi = basis(2, /*|10>*/ 0b10);
        apply_controlled_1q(Hg, psi, /*control=*/1, /*target=*/0);

        const double s = 1.0/std::sqrt(2.0);
        std::vector<C> ref(4, C{0,0});
        ref[0b10] = C{s,0}; // |10>
        ref[0b11] = C{s,0}; // |11>
        expect_state_eq(psi, ref);
    }

    // 入力 |00>（q1=0）。control=0 なので target には何も起きない
    {
        State psi = basis(2, /*|00>*/ 0b00);
        apply_controlled_1q(Hg, psi, /*control=*/1, /*target=*/0);
        State ref = basis(2, 0b00);
        expect_state_eq(psi, ref);
    }
}
