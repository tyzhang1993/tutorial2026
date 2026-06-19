#ifndef VIBRATIONALANALYSIS_H
#define VIBRATIONALANALYSIS_H

#include "Analysis.h"
#include <vector>

class Molecule;
class Hessian;
class MassProvider;

class VibrationalAnalysis : public Analysis {
public:
    VibrationalAnalysis(const Molecule& molecule, const Hessian& hessian, const MassProvider& masses);

    void run() override;
    void print(std::ostream& os) const override;

    const std::vector<double>& eigenvalues() const;
    const std::vector<double>& frequencies() const;

private:
    const Molecule& molecule_;
    const Hessian& hessian_;
    const MassProvider& masses_;

    std::vector<double> eigenvalues_;
    std::vector<double> frequencies_;
    std::vector<bool> imaginary_;
    bool has_run_ = false;
};

#endif
