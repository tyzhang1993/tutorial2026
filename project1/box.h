#include<iostream>
#include<vector>
#include<cmath>

struct vec{
    double x, y, z;
};

struct atom{
    vec v;
    int Z;
};

vec Cross(const vec& a, const vec& b);
double Dot(const vec& a, const vec& b);
double Mag(const vec& ij);
vec Vec(const vec& i, const vec& j);
vec uVec(const vec& i, const vec& j);
double bLength(const vec& i, const vec& j);
double bAngle(const vec& i, const vec& j, const vec& k);
double oAngle(const vec& i, const vec& j, const vec& k, const vec& l);
double tAngle(const vec& i, const vec& j, const vec& k, const vec& l);
vec cMass(const std::vector<atom>& atoms);
vec cMass(const std::vector<atom>& atoms);
vec Minertia(const std::vector<atom>& atoms);
double Amass(const int& Z);