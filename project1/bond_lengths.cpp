#include"box.h"

double bLength(const vec& i, const vec& j){
    vec ij = Vec(i, j);
    double lij = Mag(ij);
    
    return lij;
}