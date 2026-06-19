#ifndef MASSPROVIDER_H
#define MASSPROVIDER_H

class Molecule;

class MassProvider {
public:
    virtual ~MassProvider() = default;
    virtual double getMass(int atomic_number) const = 0;
};

#endif
