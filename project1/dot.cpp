#include"box.h"

double Dot(const vec& a, const vec& b){
    return {a.x * b.x + a.y * b.y + a.z * b.z};
}