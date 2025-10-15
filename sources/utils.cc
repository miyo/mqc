#include "qc.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <algorithm>

namespace qc{

static std::string bitstr(std::uint64_t x, int n_qubits) {
    std::string s(n_qubits, '0');
    for (int i=0;i<n_qubits;++i) if (x & (1ull<<i)) s[n_qubits-1-i] = '1'; // MSB左
    return s;
}

static std::string fmt_complex(const C& z, int prec=6) {
    std::ostringstream oss;
    oss.setf(std::ios::fixed); oss << std::setprecision(prec);
    oss << "(" << z.real();
    double im = z.imag();
    if (im >= 0) oss << "+";
    oss << im << "i)";
    return oss.str();
}

static double phase_arg(const C& z) {
    return std::atan2(z.imag(), z.real()); // [-pi, pi]
}

// 状態 |ψ> を読みやすく出力。確率降順に並べ、しきい値以下は省略。
// max_terms>0 なら上位 max_terms のみ表示。
// show_prob: 確率も表示、show_phase: 位相（ラジアン）も表示。
void pretty_print(const State& psi,
                  int n_qubits,
                  int max_terms = 0,
                  double cutoff = 1e-12,
                  int precision = 6,
                  bool show_prob = true,
                  bool show_phase = false)
{
    // 規格化チェック（軽く補正）
    double s = 0.0;
    for (auto& a : psi) s += std::norm(a);
    if (s == 0.0) {
        std::cout << "|ψ> = (all zero)\n";
        return;
    }
    s = std::sqrt(s);

    struct Item { std::uint64_t idx; C amp; double prob; };
    std::vector<Item> items; items.reserve(psi.size());
    for (std::uint64_t i=0; i<psi.size(); ++i) {
        C a = psi[i] / s;               // 軽く正規化
        double p = std::norm(a);
        if (p >= cutoff) items.push_back({i, a, p});
    }

    // 確率降順
    std::sort(items.begin(), items.end(),
              [](const Item& x, const Item& y){ return x.prob > y.prob; });

    if (max_terms > 0 && (int)items.size() > max_terms)
        items.resize(max_terms);

    // 見出し
    std::cout << "|ψ> (n=" << n_qubits << " qubits)  nonzero terms: "
              << items.size()
              << "  (cutoff=" << cutoff << ")\n";

    // 本体
    std::cout.setf(std::ios::fixed); std::cout << std::setprecision(precision);
    for (auto& it : items) {
        std::string ket = bitstr(it.idx, n_qubits);
        std::cout << "  |" << ket << ">  "
                  << "amp=" << fmt_complex(it.amp, precision);
        if (show_prob)  std::cout << "  P=" << it.prob;
        if (show_phase) std::cout << "  phase=" << phase_arg(it.amp);
        std::cout << "\n";
    }
}

}
