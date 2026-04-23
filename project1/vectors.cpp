#include"box.h"

vec Vec(const vec& i, const vec& j){
    vec v;
    v.x = j.x - i.x;
    v.y = j.y - i.y;
    v.z = j.z - i.z;
    
    return v;
}