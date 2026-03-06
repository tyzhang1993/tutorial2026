#include "VibrationalAnalysis.h"

#include "Hessian.h"
#include "MassProvider.h"
#include "Molecule.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <iomanip>
#include <stdexcept>

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace {
constexpr double kConv = 5140.484;
constexpr double kTol = 1.0e-15;
}

VibrationalAnalysis::VibrationalAnalysis(const Molecule& molecule,
                                         const Hessian& hessian,
                                         const MassProvider& masses)
    : molecule_(molecule), hessian_(hessian), masses_(masses) {}

void VibrationalAnalysis::run() {
    const int natom = molecule_.natom();
    const int ncoord = hessian_.ncoord();
    if (ncoord != 3 * natom) {
        throw std::runtime_error("Geometry/Hessian size mismatch in VibrationalAnalysis::run");
    }

    MatrixXd h_mw(ncoord, ncoord);
    for (int i = 0; i < ncoord; ++i) {
        const int atom_i = i / 3;
        const double mi = masses_.getMass(molecule_.atomicNumber(atom_i));

        for (int j = 0; j < ncoord; ++j) {
            const int atom_j = j / 3;
            const double mj = masses_.getMass(molecule_.atomicNumber(atom_j));
            h_mw(i, j) = hessian_.at(i, j) / std::sqrt(mi * mj);
        }
    }

    MatrixXd h_sym = h_mw;
    for (int i = 0; i < ncoord; ++i) {
        for (int j = i + 1; j < ncoord; ++j) {
            const double avg = 0.5 * (h_sym(i, j) + h_sym(j, i));
            h_sym(i, j) = avg;
            h_sym(j, i) = avg;
        }
    }

    Eigen::SelfAdjointEigenSolver<MatrixXd> solver(h_sym);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Eigenvalue decomposition failed");
    }

    const VectorXd evals = solver.eigenvalues();

    eigenvalues_.assign(ncoord, 0.0);
    frequencies_.assign(ncoord, 0.0);
    imaginary_.assign(ncoord, false);

    for (int k = 0; k < ncoord; ++k) {
        const double lam = evals[k];
        eigenvalues_[k] = lam;

        if (lam > kTol) {
            frequencies_[k] = std::sqrt(lam) * kConv;
            imaginary_[k] = false;
        } else if (lam < -kTol) {
            frequencies_[k] = std::sqrt(-lam) * kConv;
            imaginary_[k] = true;
        } else {
            frequencies_[k] = 0.0;
            imaginary_[k] = false;
        }
    }

    has_run_ = true;
}

void VibrationalAnalysis::print(std::ostream& os) const {
    if (!has_run_) {
        throw std::runtime_error("Call run() before print() in VibrationalAnalysis");
    }

    const int ncoord = static_cast<int>(eigenvalues_.size());

    os << "Hessian eigenvalues (hartree/amu-bohr^2):\n";
    os << std::fixed << std::setprecision(10);
    for (int k = ncoord - 1; k >= 0; --k) {
        os << std::setw(4) << k << std::setw(22) << eigenvalues_[k] << '\n';
    }

    os << '\n';
    os << "Harmonic vibrational frequencies (cm^-1):\n";
    os << std::fixed << std::setprecision(4);
    for (int k = ncoord - 1; k >= 0; --k) {
        os << std::setw(4) << k << std::setw(10) << frequencies_[k];
        if (imaginary_[k]) {
            os << 'i';
        }
        os << '\n';
    }
}

const std::vector<double>& VibrationalAnalysis::eigenvalues() const { return eigenvalues_; }
const std::vector<double>& VibrationalAnalysis::frequencies() const { return frequencies_; }
