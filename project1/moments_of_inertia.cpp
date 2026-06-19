#include "box.h"
#include <iostream>

vec Eigenv(double A[3][3]);   

vec Minertia(const std::vector<atom>& atoms){
    double Ixx = 0.0, Iyy = 0.0, Izz = 0.0;
    double Ixy = 0.0, Ixz = 0.0, Iyz = 0.0;

    vec mcen = cMass(atoms);

    for(int i = 0; i < (int)atoms.size(); i++){
        double x = atoms[i].v.x - mcen.x;
        double y = atoms[i].v.y - mcen.y;
        double z = atoms[i].v.z - mcen.z;

        double m = Amass(atoms[i].Z);

        Ixx += m * (y*y + z*z);
        Iyy += m * (x*x + z*z);
        Izz += m * (x*x + y*y);

        Ixy -= m * x * y;
        Ixz -= m * x * z;   
        Iyz -= m * y * z;  
    }

    double I[3][3] = {
        {Ixx, Ixy, Ixz},
        {Ixy, Iyy, Iyz},
        {Ixz, Iyz, Izz}
    };

    vec evI = Eigenv(I);     
    return evI;
}