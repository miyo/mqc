#include "qc.h"

namespace qc{

// ========== 例：基本ゲート ==========
void gate_X(C U[2][2]) {
    U[0][0]=C{0,0}; U[0][1]=C{1,0};
    U[1][0]=C{1,0}; U[1][1]=C{0,0};
}
void gate_H(C U[2][2]) {
    const double s = 1.0/std::sqrt(2);
    U[0][0]=C{s,0}; U[0][1]=C{s,0};
    U[1][0]=C{s,0}; U[1][1]=C{-s,0};
}
// Rz(θ) = diag(e^{-iθ/2}, e^{+iθ/2})
void gate_Rz(C U[2][2], double theta) {
    U[0][1]=U[1][0]=C{0,0};
    U[0][0]=std::exp(C{0,-theta/2});
    U[1][1]=std::exp(C{0, theta/2});
}
// CNOT 行列（行優先）
void gate_CNOT(C U4[4][4]) {
    for(int i=0;i<4;++i)for(int j=0;j<4;++j) U4[i][j]=C{0,0};
    U4[0][0]=U4[1][1]=U4[2][3]=U4[3][2]=C{1,0};
}

}
