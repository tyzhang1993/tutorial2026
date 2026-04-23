#ifndef MASS_H
#define MASS_H

#include "MassProvider.h"
#include <vector>

class Molecule;

class Mass : public MassProvider {
public:
    Mass();

    double getMass(int atomic_number) const override;
    double getMass(const Molecule& molecule, int atom_index) const;

private:
    std::vector<double> mass_table_;
};

#endif
