#include"box.h"

double Amass(const int& Z){
    std::vector<double> atomicMass = {
    0.0,      // index = Z
    1.008,    // 1  H
    4.0026,   // 2  He
    6.94,     // 3  Li
    9.0122,   // 4  Be
    10.81,    // 5  B
    12.011,   // 6  C
    14.007,   // 7  N
    15.999,   // 8  O
    18.998,   // 9  F
    20.180,   // 10 Ne
    22.990,   // 11 Na
    24.305,   // 12 Mg
    26.982,   // 13 Al
    28.085,   // 14 Si
    30.974,   // 15 P
    32.06,    // 16 S
    35.45,    // 17 Cl
    39.948,   // 18 Ar
    39.098,   // 19 K
    40.078    // 20 Ca
    };
    int m = atomicMass[Z];
    return m;
}