#include"box.h"

vec Cross(const vec& a, const vec& b){
    vec c = {
        a.y * b.z - a.z * b.y,
        a.z * b.x - b.z * a.x,
        a.x * b.y - a.y * b.x
    };
    
    return c;
}