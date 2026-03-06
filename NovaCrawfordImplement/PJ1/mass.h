#ifndef MASS_H
#define MASS_H

#include <vector>

class Mass
{
public:
    std::vector<double> mass;

    Mass();
    double printmass(int i) const;
};

#endif