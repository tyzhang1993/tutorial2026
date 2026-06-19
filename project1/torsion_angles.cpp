#include"box.h"

double tAngle(const vec& i, const vec& j, const vec& k, const vec& l){
    vec ucross_kj_kl = Cross(uVec(k, j),uVec(k, l));
    vec ucross_ji_jk = Cross(uVec(j, i),uVec(j, k));
    
    double sin_phi_jkl;
    vec kj = Vec(k, j);
    vec kl = Vec(k, l);
    double mkj = Mag(kj);
    double mkl = Mag(kl);
    double mcjkl = Mag(Cross(kj, kl));
    sin_phi_jkl = mcjkl / (mkj * mkl);
    
    double sin_phi_ijk;
    vec ji = Vec(j, i);
    vec jk = Vec(j, k);
    double mji = Mag(ji);
    double mjk = Mag(jk);
    double mcijk = Mag(Cross(ji, jk));
    sin_phi_ijk = mcijk / (mji * mjk);
    
    double cos_tau_ijkl;
    cos_tau_ijkl = Dot(ucross_ji_jk, ucross_kj_kl) / (sin_phi_ijk * sin_phi_jkl);
    
    return cos_tau_ijkl;
    
}