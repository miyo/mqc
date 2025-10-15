// sources/main_surface.cc
#include "qc.h"
#include "surface_code.h"
#include <iostream>
#include <vector>
#include <random>
#include <cstring>
#include <cstdlib>
#include <numbers>

using namespace qc;
using namespace qc::surface;

namespace {

// ---------- CLI utils ----------
bool parse_next_int(int argc, char** argv, int& i, int& out) {
    if (i + 1 >= argc) return false;
    char* endp = nullptr;
    long v = std::strtol(argv[i + 1], &endp, 10);
    if (endp == argv[i + 1] || *endp != '\0') return false;
    out = static_cast<int>(v);
    i += 1;
    return true;
}
bool parse_next_double(int argc, char** argv, int& i, double& out) {
    if (i + 1 >= argc) return false;
    char* endp = nullptr;
    double v = std::strtod(argv[i + 1], &endp);
    if (endp == argv[i + 1] || *endp != '\0') return false;
    out = v;
    i += 1;
    return true;
}
void usage(const char* prog) {
    std::cerr <<
        "Usage: " << prog << " [options]\n"
        "  --x <i>        inject X on data qubit i (0..8). Can repeat.\n"
        "  --z <i>        inject Z on data qubit i (0..8). Can repeat.\n"
        "  --y <i>        inject Y on data qubit i (0..8). Can repeat.\n"
        "  --rounds <N>   run N rounds (default: 1).\n"
        "  --noise-p <p>  depolarizing per data qubit with prob p (X/Y/Z equally).\n"
        "  --seed <u64>   RNG seed (default: random_device).\n"
        "  --help         show this help.\n";
}
inline void check_data_range(int q) {
    if (q < 0 || q > 8) {
        std::cerr << "Error: data qubit index must be in 0..8 (got " << q << ")\n";
        std::exit(2);
    }
}

// ---------- noise injection ----------
struct PauliGates {
    C Xg[2][2];
    C Zg[2][2];
    PauliGates() {
        gate_X(Xg);
        gate_Rz(Zg, std::numbers::pi); // Z (global phase ignored)
    }
};
inline void apply_pauli(int kind /*0:X,1:Z,2:Y*/, State& psi, int q, const PauliGates& G) {
    switch (kind) {
        case 0: apply_1q(G.Xg, psi, q); break;                       // X
        case 1: apply_1q(G.Zg, psi, q); break;                       // Z
        case 2: apply_1q(G.Xg, psi, q); apply_1q(G.Zg, psi, q); break; // Y = i XZ
    }
}

// 固定注入（CLI指定）＋デポラ化ノイズを適用
void inject_fixed_and_noise(State& psi,
                            const std::vector<int>& xs,
                            const std::vector<int>& zs,
                            const std::vector<int>& ys,
                            double p_noise,
                            std::mt19937_64& rng,
                            const PauliGates& G)
{
    // 固定注入
    for (int q : xs) { check_data_range(q); apply_1q(G.Xg, psi, q); }
    for (int q : zs) { check_data_range(q); apply_1q(G.Zg, psi, q); }
    for (int q : ys) { check_data_range(q); apply_1q(G.Xg, psi, q); apply_1q(G.Zg, psi, q); }

    // ノイズ（各データに独立に確率pで発生、X/Z/Y均等）
    if (p_noise > 0.0) {
        std::bernoulli_distribution coin(p_noise);
        std::uniform_int_distribution<int> which(0, 2); // 0:X,1:Z,2:Y
        for (int q = 0; q < 9; ++q) {
            if (coin(rng)) {
                int k = which(rng);
                apply_pauli(k, psi, q, G);
            }
        }
    }
}

} // namespace

int main(int argc, char** argv) {
    std::vector<int> xs, zs, ys;
    int rounds = 1;
    double p_noise = 0.0;
    bool have_seed = false;
    std::uint64_t seed = 0;

    // --- parse CLI ---
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--help") == 0) {
            usage(argv[0]);
            return 0;
        } else if (std::strcmp(argv[i], "--x") == 0) {
            int q; if (!parse_next_int(argc, argv, i, q)) { usage(argv[0]); return 1; }
            xs.push_back(q);
        } else if (std::strcmp(argv[i], "--z") == 0) {
            int q; if (!parse_next_int(argc, argv, i, q)) { usage(argv[0]); return 1; }
            zs.push_back(q);
        } else if (std::strcmp(argv[i], "--y") == 0) {
            int q; if (!parse_next_int(argc, argv, i, q)) { usage(argv[0]); return 1; }
            ys.push_back(q);
        } else if (std::strcmp(argv[i], "--rounds") == 0) {
            if (!parse_next_int(argc, argv, i, rounds) || rounds <= 0) {
                std::cerr << "Error: --rounds must be positive integer\n";
                return 1;
            }
        } else if (std::strcmp(argv[i], "--noise-p") == 0) {
            if (!parse_next_double(argc, argv, i, p_noise) || p_noise < 0.0 || p_noise > 1.0) {
                std::cerr << "Error: --noise-p must be in [0,1]\n";
                return 1;
            }
        } else if (std::strcmp(argv[i], "--seed") == 0) {
            int tmp; if (!parse_next_int(argc, argv, i, tmp)) { usage(argv[0]); return 1; }
            have_seed = true;
            seed = static_cast<std::uint64_t>(static_cast<long long>(tmp));
        } else {
            std::cerr << "Unknown option: " << argv[i] << "\n";
            usage(argv[0]);
            return 1;
        }
    }

    // RNG
    std::mt19937_64 rng(have_seed ? seed : std::random_device{}());
    PauliGates G;

    // Print header
    if (rounds > 1 || p_noise > 0.0) {
        std::cout << "# rounds=" << rounds << " noise_p=" << p_noise;
        if (have_seed) std::cout << " seed=" << seed;
        std::cout << "\n";
    }

    for (int r = 1; r <= rounds; ++r) {
        // ---- Independent run for Z syndrome ----
        State psiZ = basis(/*n=*/13, /*index=*/0);
        inject_fixed_and_noise(psiZ, xs, zs, ys, p_noise, rng, G);
        auto z = z_round(psiZ);

        // ---- Independent run for X syndrome ----
        State psiX = basis(/*n=*/13, /*index=*/0);
        prepare_all_plus_unitary(psiX); // 決定的にするために先に |+>^9
        inject_fixed_and_noise(psiX, xs, zs, ys, p_noise, rng, G);
        auto x = x_round(psiX);

        // output: 1行/ラウンド
        if (rounds == 1 && p_noise == 0.0) {
            // 旧表示（互換）
            std::cout << "Z syndrome: " << z[0] << " " << z[1] << "\n";
            std::cout << "X syndrome: " << x[0] << " " << x[1] << "\n";
        } else {
            std::cout << "round " << r << ": Z " << z[0] << " " << z[1]
                      << " | X " << x[0] << " " << x[1] << "\n";
        }
    }

    return 0;
}
