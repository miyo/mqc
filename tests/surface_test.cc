#include "qc.h"
#include "surface_code.h"
#include <gtest/gtest.h>
#include <array>
#include <numbers>
#include <iostream>

using namespace qc;
using namespace qc::surface;

// test-only: print syndromes
static void dump_syn(const char* tag, const std::array<int,2>& z, const std::array<int,2>& x) {
    std::cerr << tag << "  Z=[" << z[0] << "," << z[1] << "]  X=[" << x[0] << "," << x[1] << "]\n";
}

TEST(SurfaceD3, NoError_AllZero) {
    State psi = basis(13, 0);            // |0>^9 on data, ancilla 0
    auto z = z_round(psi);
    EXPECT_EQ(z, (std::array<int,2>{0,0}));

    prepare_all_plus_unitary(psi);
    auto x = x_round(psi);               // data -> |+>^9 にしてから X 測定
    dump_syn("NoErr", z, x);
    EXPECT_EQ(x, (std::array<int,2>{0,0}));
}

TEST(SurfaceD3, XErrorAtCenter_Z11_X00) {
    C Xg[2][2]; gate_X(Xg);
    // --- Z syndrome (independent run) ---
    State psiZ = basis(13, 0);
    apply_1q(Xg, psiZ, 4);                  // inject X
    auto z = z_round(psiZ);
    EXPECT_EQ(z, (std::array<int,2>{1,1}));
    
    // --- X syndrome (independent run) ---
    State psiX = basis(13, 0);
    prepare_all_plus_unitary(psiX);                  // H before injecting X (avoid H X H = Z)
    apply_1q(Xg, psiX, 4);                   // inject X
    auto x = x_round(psiX);
    dump_syn("X@center", z, x);
    EXPECT_EQ(x, (std::array<int,2>{0,0}));
}

TEST(SurfaceD3, ZErrorAtCenter_X11_Z00) {
    C Zg[2][2]; gate_Rz(Zg, std::numbers::pi);
    // --- X syndrome (independent run) ---
    State psiX = basis(13, 0);
    prepare_all_plus_unitary(psiX);                  // H before injecting Z (H Z H = X にはならない)
    apply_1q(Zg, psiX, 4);                   // inject Z
    auto x = x_round(psiX);                  // expect [1,1]
    dump_syn("Z@center(Xrun)", std::array<int,2>{0,0}, x);
    EXPECT_EQ(x, (std::array<int,2>{1,1}));
    
    // --- Z syndrome (independent run) ---
    State psiZ = basis(13, 0);
    apply_1q(Zg, psiZ, 4);
    auto z = z_round(psiZ);                  // expect [0,0]
    EXPECT_EQ(z, (std::array<int,2>{0,0}));
}

TEST(SurfaceD3, YErrorAtCenter_Both11) {
    C Xg[2][2]; gate_X(Xg);
    C Zg[2][2]; gate_Rz(Zg, std::numbers::pi);
    // --- Z syndrome (independent run) ---
    State psiZ = basis(13, 0);
    apply_1q(Xg, psiZ, 4);                   // Y の X 成分だけで Z は反応
    apply_1q(Zg, psiZ, 4);                   // 形式的に Y 完成
    auto z = z_round(psiZ);                  // expect [1,1]
    EXPECT_EQ(z, (std::array<int,2>{1,1}));
    
    // --- X syndrome (independent run) ---
    State psiX = basis(13, 0);
    prepare_all_plus_unitary(psiX);
    apply_1q(Xg, psiX, 4);
    apply_1q(Zg, psiX, 4);
    auto x = x_round(psiX);                  // expect [1,1]
    dump_syn("Y@center", z, x);
}
