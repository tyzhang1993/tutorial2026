#include "mass.h"

Mass::Mass()
{
    mass.push_back(0.0000000);
    mass.push_back(1.007825);
    mass.push_back(4.002602);
    mass.push_back(7.01600);
    mass.push_back(9.012182);
    mass.push_back(11.00931);
    mass.push_back(12.00000);
    mass.push_back(14.00307);
    mass.push_back(15.99491);
    mass.push_back(18.99840);
    mass.push_back(19.99244);
}

double Mass::printmass(int i) const
{
    return mass[i];
}