#pragma once

#include <vector>
#include <complex>
#include <cstdint>
#include <cmath>

namespace qc{

using C     = std::complex<double>;
using State = std::vector<C>;

State basis(int n_qubits, std::uint64_t index);
void renormalize(State& psi);

std::uint64_t measure_all(State& psi);
int measure_qubit_Z(State& psi, int target);

void apply_1q(const C U[2][2], State& psi, int target);
void apply_2q(const C U4[4][4], State& psi, int qA, int qB);
void apply_controlled_1q(const C U[2][2], State& psi, int control, int target);

void gate_X(C U[2][2]);
void gate_H(C U[2][2]);
void gate_Rz(C U[2][2], double theta);
void gate_CNOT(C U4[4][4]);

void pretty_print(const State& psi, int n_qubits, int max_terms, double cutoff, int precision, bool show_prob, bool show_phase);

}
