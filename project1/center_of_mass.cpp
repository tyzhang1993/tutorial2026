#include"box.h"

vec cMass(const std::vector<atom>& atoms){
    double sum_m = 0.0, 
           sum_mx = 0.0, 
           sum_my = 0.0, 
           sum_mz = 0.0;
    
    for(int i = 0; i < atoms.size(); i++){
        double m = Amass(atoms[i].Z);
        sum_m += m;
        sum_mx += m * atoms[i].v.x;
        sum_my += m * atoms[i].v.y;
        sum_mz += m * atoms[i].v.z;
    }
    
    vec cmass = {sum_mx/sum_m, sum_my/sum_m, sum_mz/sum_m};
    return cmass;
}