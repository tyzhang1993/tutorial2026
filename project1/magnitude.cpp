#include"box.h"

double Mag(const vec& ij){
    return {std::sqrt(Dot(ij, ij))};
}