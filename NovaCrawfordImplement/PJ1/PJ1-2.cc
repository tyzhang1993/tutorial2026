#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

#include "molecule.h"
#include "mass.h"

using namespace std;

static double det3(double A[3][3])
{
    return
        A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
        A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
        A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
}

static void eigenvalues_symm_3x3(double A[3][3], double eval[3])
{
    double p1 = A[0][1] * A[0][1] + A[0][2] * A[0][2] + A[1][2] * A[1][2];

    if (p1 < 1.0e-16)
    {
        eval[0] = A[0][0];
        eval[1] = A[1][1];
        eval[2] = A[2][2];
    }
    else
    {
        double q = (A[0][0] + A[1][1] + A[2][2]) / 3.0;

        double B[3][3];
        double p2 =
            (A[0][0] - q) * (A[0][0] - q) +
            (A[1][1] - q) * (A[1][1] - q) +
            (A[2][2] - q) * (A[2][2] - q) +
            2.0 * p1;

        double p = std::sqrt(p2 / 6.0);

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                B[i][j] = (A[i][j] - (i == j ? q : 0.0)) / p;

        double r = det3(B) / 2.0;
        if (r > 1.0) r = 1.0;
        if (r < -1.0) r = -1.0;

        double phi = std::acos(r) / 3.0;

        eval[0] = q + 2.0 * p * std::cos(phi);
        eval[2] = q + 2.0 * p * std::cos(phi + 2.0 * M_PI / 3.0);
        eval[1] = 3.0 * q - eval[0] - eval[2];
    }

    if (eval[0] > eval[1]) std::swap(eval[0], eval[1]);
    if (eval[1] > eval[2]) std::swap(eval[1], eval[2]);
    if (eval[0] > eval[1]) std::swap(eval[0], eval[1]);
}

int main()
{
    Molecule mol("geom.dat", 0);
    Mass m;

    cout << "Number of atoms: " << mol.natom << endl;
    cout << "Input Cartesian coordinates:" << endl;
    mol.print_geom();

    cout << "\nInteratomic distances (bohr):" << endl;
    for (int i = 0; i < mol.natom; i++)
        for (int j = 0; j < i; j++)
            cout << i << " " << j << "  "
                 << fixed << setprecision(5)
                 << mol.bond(i, j) << endl;

    cout << "\nBond angles:\n" << endl;
    for (int i = 0; i < mol.natom; i++)
        for (int j = 0; j < i; j++)
            for (int k = 0; k < j; k++)
                if (mol.bond(i, j) < 4.0 && mol.bond(j, k) < 4.0)
                    cout << setw(2) << k << "-" << setw(2) << j << "-" << setw(2) << i << " "
                         << fixed << setprecision(6)
                         << mol.angle(k, j, i) * 180.0 / M_PI << endl;

    cout << "\nOut-of-plane angles:\n" << endl;
    for (int i = 0; i < mol.natom; i++)
        for (int j = 0; j < mol.natom; j++)
            for (int k = 0; k < mol.natom; k++)
                for (int l = 0; l < k; l++)
                    if (i != j && i != k && i != l &&
                        j != k && j != l && k != l &&
                        mol.bond(i, j) < 4.0 &&
                        mol.bond(j, k) < 4.0 &&
                        mol.bond(k, l) < 4.0)
                        cout << setw(2) << i << "-" << setw(2) << j << "-" << setw(2) << k << "-" << setw(2) << l << " "
                             << fixed << setprecision(6)
                             << mol.oop(i, k, j, l) * 180.0 / M_PI << endl;

    cout << "\nTorsional angles:\n" << endl;
    for (int i = 0; i < mol.natom; i++)
        for (int j = 0; j < i; j++)
            for (int k = 0; k < j; k++)
                for (int l = 0; l < k; l++)
                    if (mol.bond(i, j) < 4.0 &&
                        mol.bond(j, k) < 4.0 &&
                        mol.bond(k, l) < 4.0)
                        cout << setw(2) << i << "-" << setw(2) << j << "-" << setw(2) << k << "-" << setw(2) << l << " "
                             << fixed << setprecision(6)
                             << mol.torsion(i, j, k, l) * 180.0 / M_PI << endl;

    double mtot = 0.0, xm = 0.0, ym = 0.0, zm = 0.0;

    for (int i = 0; i < mol.natom; i++)
        mtot += m.printmass(mol.zvals[i]);

    for (int i = 0; i < mol.natom; i++)
    {
        double mi = m.printmass(mol.zvals[i]);
        xm += mi * mol.geom[i][0];
        ym += mi * mol.geom[i][1];
        zm += mi * mol.geom[i][2];
    }

    xm /= mtot;
    ym /= mtot;
    zm /= mtot;

    cout << "\nMolecular center of mass: "
         << setw(12) << fixed << setprecision(8) << xm
         << setw(12) << ym
         << setw(12) << zm << endl;

    mol.translate(-xm, -ym, -zm);

    double I[3][3] = {{0.0, 0.0, 0.0},
                      {0.0, 0.0, 0.0},
                      {0.0, 0.0, 0.0}};

    for (int i = 0; i < mol.natom; i++)
    {
        double mi = m.printmass(mol.zvals[i]);
        double x = mol.geom[i][0];
        double y = mol.geom[i][1];
        double z = mol.geom[i][2];

        I[0][0] += mi * (y * y + z * z);
        I[1][1] += mi * (x * x + z * z);
        I[2][2] += mi * (x * x + y * y);

        I[0][1] -= mi * x * y;
        I[0][2] -= mi * x * z;
        I[1][2] -= mi * y * z;
    }

    I[1][0] = I[0][1];
    I[2][0] = I[0][2];
    I[2][1] = I[1][2];

    cout << "\nMoment of inertia tensor:\n" << endl;
    cout << "           1           2           3\n" << endl;
    for (int r = 0; r < 3; r++)
    {
        cout << "    " << r + 1;
        for (int c = 0; c < 3; c++)
            cout << setw(12) << fixed << setprecision(7) << I[r][c];
        cout << endl;
    }

    double eval[3];
    eigenvalues_symm_3x3(I, eval);

    cout << "\nPrincipal moments of inertia (amu * bohr^2):\n";
    cout << "\t" << fixed << setprecision(6)
         << setw(10) << eval[0] << "  "
         << setw(10) << eval[1] << "  "
         << setw(10) << eval[2] << endl;

    const double bohr_to_A = 0.529177249;
    cout << "\nPrincipal moments of inertia (amu * AA^2):\n";
    cout << "\t" << fixed << setprecision(6)
         << setw(10) << eval[0] * bohr_to_A * bohr_to_A << "  "
         << setw(10) << eval[1] * bohr_to_A * bohr_to_A << "  "
         << setw(10) << eval[2] * bohr_to_A * bohr_to_A << endl;

    const double amu_to_g = 1.6605402e-24;
    const double bohr_to_cm = bohr_to_A * 1.0e-8;
    cout << "\nPrincipal moments of inertia (g * cm^2):\n";
    cout << "\t" << scientific << setprecision(6)
         << eval[0] * amu_to_g * bohr_to_cm * bohr_to_cm << " "
         << eval[1] * amu_to_g * bohr_to_cm * bohr_to_cm << " "
         << eval[2] * amu_to_g * bohr_to_cm * bohr_to_cm << endl;

    cout << "\nMolecule is ";
    if (mol.natom == 2) cout << "a diatomic molecule." << endl;
    else if (fabs(eval[0]) < 1.0e-4) cout << "linear." << endl;
    else if (fabs(eval[0] - eval[1]) < 1.0e-4 && fabs(eval[1] - eval[2]) < 1.0e-4) cout << "a spherical top." << endl;
    else if (fabs(eval[0] - eval[1]) < 1.0e-4 && fabs(eval[1] - eval[2]) > 1.0e-4) cout << "an oblate symmetric top." << endl;
    else if (fabs(eval[0] - eval[1]) > 1.0e-4 && fabs(eval[1] - eval[2]) < 1.0e-4) cout << "a prolate symmetric top." << endl;
    else cout << "an asymmetric top." << endl;

    const double MHz_factor = 1804741.132488276;
    double A = MHz_factor / eval[0];
    double B = MHz_factor / eval[1];
    double C = MHz_factor / eval[2];

    cout << "\nRotational constants (MHz):" << endl;
    cout << "\tA = " << fixed << setprecision(3) << A
         << "\t B = " << B
         << "\t C = " << C << endl;

    const double MHz_per_cm1 = 29979.2458;
    cout << "\nRotational constants (cm-1):" << endl;
    cout << "\tA = " << fixed << setprecision(4) << A / MHz_per_cm1
         << "\t B = " << B / MHz_per_cm1
         << "\t C = " << C / MHz_per_cm1 << endl;

    return 0;
}