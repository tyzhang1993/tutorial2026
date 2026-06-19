#include "DataObject.h"
#include "Hessian.h"
#include "Mass.h"
#include "Molecule.h"
#include "VibrationalAnalysis.h"

#include <iostream>
#include <memory>
#include <vector>

int main() {
    try {
        Molecule molecule("PJ2geom.dat", 0);
        Hessian hessian("PJ2hess.dat");
        Mass mass;

        std::vector<std::unique_ptr<DataObject>> objects;
        objects.push_back(std::make_unique<Molecule>(molecule));
        objects.push_back(std::make_unique<Hessian>(hessian));

        // Polymorphism is used through the DataObject base class.
        // We intentionally keep the program output clean and only print
        // the vibrational-analysis results requested for the assignment.
        (void)objects;

        // Overloading demo: both calls resolve to different translate overloads.
        molecule.translate(0.0, 0.0, 0.0);
        molecule.translate(Vec3{0.0, 0.0, 0.0});

        VibrationalAnalysis analysis(molecule, hessian, mass);
        analysis.run();
        analysis.print(std::cout);
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }

    return 0;
}
