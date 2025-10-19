#include "qc.h"
#include "surface_code.h"
#include <gtest/gtest.h>
#include <array>
#include <numbers>
#include <iostream>

using namespace qc;
using namespace qc::surface;

// test-only: print syndromes
static void dump_syn(const char* tag, const std::vector<int>& z, const std::vector<int>& x) {
    std::cerr << tag << "  Z=["; 
    for (size_t i=0;i<z.size();++i){ if(i)std::cerr<<","; std::cerr<<z[i]; }
    std::cerr << "]  X=[";
    for (size_t i=0;i<x.size();++i){ if(i)std::cerr<<","; std::cerr<<x[i]; }
    std::cerr << "]\n";
}

TEST(SurfaceD3, NoError_AllZero) {
    auto sc = build_surface_code(3);
    State psi = basis(sc.n_qubits(), 0);            // |0>^9 on data, ancilla 0
    auto z = z_round(psi, sc);
    EXPECT_EQ(z, (std::vector<int>({0,0})));

    prepare_all_plus_unitary(psi, sc);
    auto x = x_round(psi, sc);               // data -> |+>^9 にしてから X 測定
    dump_syn("NoErr", z, x);
    EXPECT_EQ(x, (std::vector<int>({0,0})));
}

TEST(SurfaceD3, XErrorAtCenter_Z11_X00) {
    auto sc = build_surface_code(3);
    C Xg[2][2]; gate_X(Xg);
    // --- Z syndrome (independent run) ---
    State psiZ = basis(sc.n_qubits(), 0);
    apply_1q(Xg, psiZ, 4);                  // inject X
    auto z = z_round(psiZ, sc);
    EXPECT_EQ(z, (std::vector<int>({1,1})));
    
    // --- X syndrome (independent run) ---
    State psiX = basis(sc.n_qubits(), 0);
    prepare_all_plus_unitary(psiX, sc);                  // H before injecting Z (note: H Z H = X)
    apply_1q(Xg, psiX, 4);                   // inject X
    auto x = x_round(psiX, sc);
    dump_syn("X@center", z, x);
    EXPECT_EQ(x, (std::vector<int>({0,0})));
}

TEST(SurfaceD3, ZErrorAtCenter_X11_Z00) {
    auto sc = build_surface_code(3);
    C Zg[2][2]; gate_Rz(Zg, std::numbers::pi);
    // --- X syndrome (independent run) ---
    State psiX = basis(sc.n_qubits(), 0);
    prepare_all_plus_unitary(psiX, sc);                  // H before injecting Z (note: H Z H = X)
    apply_1q(Zg, psiX, 4);                   // inject Z
    auto x = x_round(psiX, sc);                  // expect [1,1]
    dump_syn("Z@center(Xrun)", std::vector<int>({0,0}), x);
    EXPECT_EQ(x, (std::vector<int>({1,1})));
    
    // --- Z syndrome (independent run) ---
    State psiZ = basis(13, 0);
    apply_1q(Zg, psiZ, 4);
    auto z = z_round(psiZ, sc);                  // expect [0,0]
    EXPECT_EQ(z, (std::vector<int>({0,0})));
}

TEST(SurfaceD3, YErrorAtCenter_Both11) {
    auto sc = build_surface_code(3);
    C Xg[2][2]; gate_X(Xg);
    C Zg[2][2]; gate_Rz(Zg, std::numbers::pi);
    // --- Z syndrome (independent run) ---
    State psiZ = basis(sc.n_qubits(), 0);
    apply_1q(Xg, psiZ, 4);                   // Y の X 成分だけで Z は反応
    apply_1q(Zg, psiZ, 4);                   // 形式的に Y 完成
    auto z = z_round(psiZ, sc);                  // expect [1,1]
    EXPECT_EQ(z, (std::vector<int>({1,1})));
    
    // --- X syndrome (independent run) ---
    State psiX = basis(sc.n_qubits(), 0);
    prepare_all_plus_unitary(psiX, sc);
    apply_1q(Xg, psiX, 4);
    apply_1q(Zg, psiX, 4);
    auto x = x_round(psiX, sc);                  // expect [1,1]
    EXPECT_EQ(x, (std::vector<int>({1,1})));
    dump_syn("Y@center", z, x);
}

