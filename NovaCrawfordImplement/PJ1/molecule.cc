#include "molecule.h"
#include <cstdio>
#include <fstream>
#include <cassert>
#include <cmath>

Molecule::Molecule(const char *filename, int q)
{
    charge = q;

    std::ifstream is(filename);
    assert(is.good());

    is >> natom;

    zvals = new int[natom];
    geom = new double *[natom];
    for (int i = 0; i < natom; i++)
        geom[i] = new double[3];

    for (int i = 0; i < natom; i++)
        is >> zvals[i] >> geom[i][0] >> geom[i][1] >> geom[i][2];

    is.close();
}

Molecule::~Molecule()
{
    delete[] zvals;
    for (int i = 0; i < natom; i++)
        delete[] geom[i];
    delete[] geom;
}

void Molecule::print_geom()
{
    for (int i = 0; i < natom; i++)
        std::printf("%d %16.12f %16.12f %16.12f\n",
                    zvals[i], geom[i][0], geom[i][1], geom[i][2]);
}

void Molecule::translate(double x, double y, double z)
{
    for (int i = 0; i < natom; i++)
    {
        geom[i][0] += x;
        geom[i][1] += y;
        geom[i][2] += z;
    }
}

void Molecule::rotate(double phi)
{
    double c = std::cos(phi);
    double s = std::sin(phi);

    for (int i = 0; i < natom; i++)
    {
        double x0 = geom[i][0];
        double y0 = geom[i][1];

        geom[i][0] = c * x0 - s * y0;
        geom[i][1] = s * x0 + c * y0;
    }
}

double Molecule::bond(int i, int j)
{
    double dx = geom[i][0] - geom[j][0];
    double dy = geom[i][1] - geom[j][1];
    double dz = geom[i][2] - geom[j][2];

    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double Molecule::unit(int cart, int a, int b)
{
    return -(geom[a][cart] - geom[b][cart]) / bond(a, b);
}

double Molecule::angle(int i, int j, int k)
{
    double cosang =
        unit(0, j, i) * unit(0, j, k) +
        unit(1, j, i) * unit(1, j, k) +
        unit(2, j, i) * unit(2, j, k);

    if (cosang > 1.0) cosang = 1.0;
    if (cosang < -1.0) cosang = -1.0;

    return std::acos(cosang);
}

double Molecule::oop(int a, int b, int c, int d)
{
    double bcdx = unit(1, c, b) * unit(2, c, d) - unit(2, c, b) * unit(1, c, d);
    double bcdy = unit(2, c, b) * unit(0, c, d) - unit(0, c, b) * unit(2, c, d);
    double bcdz = unit(0, c, b) * unit(1, c, d) - unit(1, c, b) * unit(0, c, d);

    double val = bcdx * unit(0, c, a) +
                 bcdy * unit(1, c, a) +
                 bcdz * unit(2, c, a);

    double s = std::sin(angle(b, c, d));
    if (std::fabs(s) < 1.0e-14)
        return 0.0;

    val /= s;

    if (val > 1.0) val = 1.0;
    if (val < -1.0) val = -1.0;

    return std::asin(val);
}

double Molecule::torsion(int i, int j, int k, int l)
{
    double e1x = unit(1, i, j) * unit(2, j, k) - unit(2, i, j) * unit(1, j, k);
    double e1y = unit(2, i, j) * unit(0, j, k) - unit(0, i, j) * unit(2, j, k);
    double e1z = unit(0, i, j) * unit(1, j, k) - unit(1, i, j) * unit(0, j, k);

    double e2x = unit(1, j, k) * unit(2, k, l) - unit(2, j, k) * unit(1, k, l);
    double e2y = unit(2, j, k) * unit(0, k, l) - unit(0, j, k) * unit(2, k, l);
    double e2z = unit(0, j, k) * unit(1, k, l) - unit(1, j, k) * unit(0, k, l);

    double denom = std::sin(angle(i, j, k)) * std::sin(angle(j, k, l));
    if (std::fabs(denom) < 1.0e-14)
        return 0.0;

    double cos_tor = (e1x * e2x + e1y * e2y + e1z * e2z) / denom;

    if (cos_tor > 1.0) cos_tor = 1.0;
    if (cos_tor < -1.0) cos_tor = -1.0;

    double tau = std::acos(cos_tor);

    double cx = e1y * e2z - e1z * e2y;
    double cy = e1z * e2x - e1x * e2z;
    double cz = e1x * e2y - e1y * e2x;

    double sign = cx * unit(0, j, k) + cy * unit(1, j, k) + cz * unit(2, j, k);
    if (sign < 0.0)
        tau = -tau;

    return tau;
}
