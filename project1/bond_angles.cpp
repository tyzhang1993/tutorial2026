#include"box.h"

double bAngle(const vec& i, const vec& j, const vec& k){
    double angle_ijk;
    vec ji = uVec(j, i);
    vec jk = uVec(j, k);
    angle_ijk = Dot(ji, jk);
    
    return angle_ijk;
}