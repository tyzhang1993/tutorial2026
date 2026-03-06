#include "Molecule.h"

#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace {
constexpr double kPi = 3.14159265358979323846;
}

Molecule::Molecule(const std::string& filename, int charge)
    : filename_(filename), charge_(charge) {
    std::ifstream is(filename);
    if (!is) {
        throw std::runtime_error("Could not open geometry file: " + filename);
    }

    int natom_read = 0;
    is >> natom_read;
    if (!is || natom_read <= 0) {
        throw std::runtime_error("Invalid atom count in geometry file: " + filename);
    }

    zvals_.resize(natom_read);
    geom_.resize(natom_read);

    for (int i = 0; i < natom_read; ++i) {
        double z = 0.0;
        is >> z >> geom_[i].x >> geom_[i].y >> geom_[i].z;
        if (!is) {
            throw std::runtime_error("Failed to read atom record from geometry file: " + filename);
        }
        zvals_[i] = static_cast<int>(std::lround(z));
    }
}

int Molecule::natom() const { return static_cast<int>(zvals_.size()); }
int Molecule::charge() const { return charge_; }
int Molecule::atomicNumber(int atom) const { return zvals_.at(atom); }
const Vec3& Molecule::coord(int atom) const { return geom_.at(atom); }

void Molecule::printSummary(std::ostream& os) const {
    os << "Molecule(" << filename_ << "), natom=" << natom() << ", charge=" << charge_;
}

std::string Molecule::className() const { return "Molecule"; }

void Molecule::translate(double dx, double dy, double dz) {
    for (auto& r : geom_) {
        r.x += dx;
        r.y += dy;
        r.z += dz;
    }
}

void Molecule::translate(const Vec3& delta) {
    translate(delta.x, delta.y, delta.z);
}

double Molecule::bond(int i, int j) const {
    const Vec3& a = geom_.at(i);
    const Vec3& b = geom_.at(j);
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double Molecule::unit(int cart, int a, int b) const {
    const Vec3& ra = geom_.at(a);
    const Vec3& rb = geom_.at(b);
    const double r = bond(a, b);
    if (r == 0.0) {
        throw std::runtime_error("Zero bond length encountered in Molecule::unit");
    }

    switch (cart) {
        case 0: return -(ra.x - rb.x) / r;
        case 1: return -(ra.y - rb.y) / r;
        case 2: return -(ra.z - rb.z) / r;
        default: throw std::out_of_range("Cartesian index must be 0, 1, or 2");
    }
}

double Molecule::clampToUnit(double x) {
    if (x < -1.0) return -1.0;
    if (x > 1.0) return 1.0;
    return x;
}

double Molecule::angle(int i, int j, int k) const {
    const double ang = unit(0, j, i) * unit(0, j, k)
                     + unit(1, j, i) * unit(1, j, k)
                     + unit(2, j, i) * unit(2, j, k);
    return std::acos(clampToUnit(ang));
}

double Molecule::oop(int a, int b, int c, int d) const {
    const double bcdx = unit(1, c, b) * unit(2, c, d) - unit(2, c, b) * unit(1, c, d);
    const double bcdy = unit(2, c, b) * unit(0, c, d) - unit(0, c, b) * unit(2, c, d);
    const double bcdz = unit(0, c, b) * unit(1, c, d) - unit(1, c, b) * unit(0, c, d);

    const double exx = bcdx * unit(0, c, a);
    const double eyy = bcdy * unit(1, c, a);
    const double ezz = bcdz * unit(2, c, a);

    const double theta = (exx + eyy + ezz) / std::sin(angle(b, c, d));
    return std::asin(clampToUnit(theta));
}

double Molecule::torsion(int i, int j, int k, int l) const {
    const double jkklx = unit(1, j, k) * unit(2, k, l) - unit(2, j, k) * unit(1, k, l);
    const double jkkly = unit(2, j, k) * unit(0, k, l) - unit(0, j, k) * unit(2, k, l);
    const double jkklz = unit(0, j, k) * unit(1, k, l) - unit(1, j, k) * unit(0, k, l);

    const double ijjkx = unit(1, i, j) * unit(2, j, k) - unit(2, i, j) * unit(1, j, k);
    const double ijjky = unit(2, i, j) * unit(0, j, k) - unit(0, i, j) * unit(2, j, k);
    const double ijjkz = unit(0, i, j) * unit(1, j, k) - unit(1, i, j) * unit(0, j, k);

    const double tor = (jkklx * ijjkx + jkkly * ijjky + jkklz * ijjkz)
                     / (std::sin(angle(i, j, k)) * std::sin(angle(j, k, l)));

    double value = std::acos(clampToUnit(tor));

    double cross_x = ijjky * jkklz - ijjkz * jkkly;
    double cross_y = ijjkz * jkklx - ijjkx * jkklz;
    double cross_z = ijjkx * jkkly - ijjky * jkklx;
    const double norm2 = cross_x * cross_x + cross_y * cross_y + cross_z * cross_z;
    if (norm2 == 0.0) {
        return value;
    }

    cross_x /= norm2;
    cross_y /= norm2;
    cross_z /= norm2;

    const double dot = cross_x * unit(0, j, k)
                     + cross_y * unit(1, j, k)
                     + cross_z * unit(2, j, k);
    const double sign = (dot < 0.0) ? -1.0 : 1.0;
    return value * sign;
}
