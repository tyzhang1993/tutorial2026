#include"box.h"

double oAngle(const vec& i, const vec& j, const vec& k, const vec& l){
    double sin_theta_ijkl;
    
    vec ucross_kj_kl;
    ucross_kj_kl = Cross(uVec(k, j),uVec(k, l));
    
    double sin_phi_jkl;
    vec kj = Vec(k, j);
    vec kl = Vec(k, l);
    double mkj = Mag(kj);
    double mkl = Mag(kl);
    double mcjkl = Mag(Cross(kj, kl));
    sin_phi_jkl = mcjkl / (mkj * mkl);
    
    sin_theta_ijkl = Dot(ucross_kj_kl, uVec(k, i)) / sin_phi_jkl;
    
    return sin_theta_ijkl;
}