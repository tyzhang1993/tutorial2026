#include "Mass.h"
#include "Molecule.h"

#include <stdexcept>

Mass::Mass()
    : mass_table_{
          0.0000000, 1.007825, 4.002602, 7.01600, 9.012182, 11.00931, 12.00000,
          14.00307, 15.99491, 18.99840, 19.99244, 22.98977, 23.98504, 26.98154,
          27.97693, 30.97376, 31.97207, 35.45300, 38.96371} {}

double Mass::getMass(int atomic_number) const {
    if (atomic_number < 0 || atomic_number >= static_cast<int>(mass_table_.size())) {
        throw std::out_of_range("Atomic number out of range in Mass::getMass");
    }
    return mass_table_[static_cast<std::size_t>(atomic_number)];
}

double Mass::getMass(const Molecule& molecule, int atom_index) const {
    return getMass(molecule.atomicNumber(atom_index));
}
