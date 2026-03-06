#ifndef MOLECULE_H
#define MOLECULE_H

#include "DataObject.h"
#include "Vec3.h"
#include <string>
#include <vector>

class Molecule : public DataObject {
public:
    Molecule(const std::string& filename, int charge);

    int natom() const;
    int charge() const;
    int atomicNumber(int atom) const;
    const Vec3& coord(int atom) const;

    void printSummary(std::ostream& os) const override;
    std::string className() const override;

    void translate(double dx, double dy, double dz);
    void translate(const Vec3& delta);

    double bond(int atom1, int atom2) const;
    double angle(int atom1, int atom2, int atom3) const;
    double torsion(int atom1, int atom2, int atom3, int atom4) const;
    double oop(int atom1, int atom2, int atom3, int atom4) const;
    double unit(int cart, int atom1, int atom2) const;

private:
    std::string filename_;
    int charge_ = 0;
    std::vector<int> zvals_;
    std::vector<Vec3> geom_;

    static double clampToUnit(double x);
};

#endif
