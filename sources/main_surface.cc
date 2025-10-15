// sources/main_surface.cc
#include "qc.h"
#include "surface_code.h"
#include <iostream>
#include <vector>
#include <cstring>
#include <cstdlib>

using namespace qc;
using namespace qc::surface;

namespace {

bool parse_next_int(int argc, char** argv, int& i, int& out) {
    if (i + 1 >= argc) return false;
    char* endp = nullptr;
    long v = std::strtol(argv[i + 1], &endp, 10);
    if (endp == argv[i + 1] || *endp != '\0') return false;
    out = static_cast<int>(v);
    i += 1;
    return true;
}

void usage(const char* prog) {
    std::cerr <<
        "Usage: " << prog << " [options]\n"
        "  --x <i>    inject X on data qubit i (0..8). Can repeat.\n"
        "  --z <i>    inject Z on data qubit i (0..8). Can repeat.\n"
        "  --y <i>    inject Y on data qubit i (0..8). Can repeat.\n"
        "  --help     show this help.\n";
}

inline void check_data_range(int q) {
    if (q < 0 || q > 8) {
        std::cerr << "Error: data qubit index must be in 0..8 (got " << q << ")\n";
        std::exit(2);
    }
}

} // namespace

int main(int argc, char** argv) {
    std::vector<int> xs, zs, ys;

    // --- parse CLI ---
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--help") == 0) {
            usage(argv[0]);
            return 0;
        } else if (std::strcmp(argv[i], "--x") == 0) {
            int q; if (!parse_next_int(argc, argv, i, q)) { usage(argv[0]); return 1; }
            check_data_range(q); xs.push_back(q);
        } else if (std::strcmp(argv[i], "--z") == 0) {
            int q; if (!parse_next_int(argc, argv, i, q)) { usage(argv[0]); return 1; }
            check_data_range(q); zs.push_back(q);
        } else if (std::strcmp(argv[i], "--y") == 0) {
            int q; if (!parse_next_int(argc, argv, i, q)) { usage(argv[0]); return 1; }
            check_data_range(q); ys.push_back(q);
        } else {
            std::cerr << "Unknown option: " << argv[i] << "\n";
            usage(argv[0]);
            return 1;
        }
    }

    // Gates
    C Xg[2][2]; gate_X(Xg);
    C Zg[2][2]; gate_Rz(Zg, std::numbers::pi); // Z up to global phase

    // ========== Independent run for Z syndrome ==========
    State psiZ = basis(/*n=*/13, /*index=*/0);
    // Inject errors as requested
    for (int q : xs) apply_1q(Xg, psiZ, q);
    for (int q : zs) apply_1q(Zg, psiZ, q);
    for (int q : ys) { apply_1q(Xg, psiZ, q); apply_1q(Zg, psiZ, q); }
    // Measure Z round
    auto z = z_round(psiZ);

    // ========== Independent run for X syndrome ==========
    State psiX = basis(/*n=*/13, /*index=*/0);
    // Make X-syndrome deterministic: prepare |+>^9 first
    prepare_all_plus_unitary(psiX);
    // Inject the *same* errors AFTER making |+>^9 so H doesn't convert them
    for (int q : xs) apply_1q(Xg, psiX, q);                       // X stays X
    for (int q : zs) apply_1q(Zg, psiX, q);                       // Z stays Z (no extra H later)
    for (int q : ys) { apply_1q(Xg, psiX, q); apply_1q(Zg, psiX, q); }
    // Measure X round
    auto x = x_round(psiX);

    // --- print results ---
    std::cout << "Z syndrome: " << z[0] << " " << z[1] << "\n";
    std::cout << "X syndrome: " << x[0] << " " << x[1] << "\n";
    return 0;
}

