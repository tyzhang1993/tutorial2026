#include"box.h"

vec uVec(const vec& i, const vec& j){
    vec ij = Vec(i, j);
    
    vec uij;
    uij.x = ij.x/Mag(ij);
    uij.y = ij.y/Mag(ij);
    uij.z = ij.z/Mag(ij);
    
    return uij;
}