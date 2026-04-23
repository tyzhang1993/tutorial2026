#include<fstream>
#include<string>
#include"box.h"

int main() {
	std::cout << "Please enter your file path: " << std::endl;
	std::string path;
	std::getline(std::cin, path);
	std::ifstream read(path);

	while(true) {
		if(path.empty()) {
			return 0;
		}

		if(!read) {
			std::cout << "Can't open file. Try again: " << std::endl;
			continue;
		} else {
			std::cout << "File opened sucessfully!" << std::endl;
			break;
		}
	}

	std::vector<atom> atoms;
	atom a{};
	while(read >> a.Z >> a.v.x >> a.v.y >> a.v.z) {
		atoms.push_back(a);
	}
	
	int s = atoms.size();

	std::cout << "Your data: " << std::endl << "Index, Z, x, y, z" << std::endl;
	for(int i = 0; i < s; i++) {
		std::cout << i + 1 << ", " << atoms[i].Z 
		          << ", " << atoms[i].v.x << ", "
		          << atoms[i].v.y << ", "
		          << atoms[i].v.z << std::endl;
	}
	
	std::cout << "------Bond Lengths------" << std::endl;
	for(int i = 0; i < s; i++){
	    for(int j = i+1; j < s; j++){
	        double l = bLength(atoms[i].v, atoms[j].v);
	        std::cout << "bond length(" << i+1 << j+1 << ")= " << l;
	    }
	}
	
	std::cout << "------Bond Angles------" << std::endl;
	for(int j = 0; j < s; j++){
	    for(int i = 0; i < s; i++){
	        if(i == j)continue;
	        for(int k = i+1; k < s; k++){
	            if(k == j)continue;
	            double angle_ijk = bAngle(atoms[i].v, atoms[j].v, atoms[k].v);
	            std::cout << "cosΦ(" << i+1 << j+1 << k+1
	                      << ")= " << angle_ijk << std::endl;
	        }
	    }
	}
	
	std::cout << "------Out of Plane Angles------" << std::endl;
	for(int k = 0; k < s; k++){
	    for(int j = 0; j < s; j++){
	        if(j == k)continue;
	        for(int l = j+1; l < s; l++){
	            if(l == k)continue;
	            for(int i = 0; i < s; i++){
	                if(i == l || i == j || i == k)continue;
	                double sin_theta_ijkl = oAngle(atoms[i].v, atoms[j].v, atoms[k].v, atoms[l].v);
	                std::cout << "sinθ(" << i+1 << j+1 << k+1 << l+1
	                      << ")= " << sin_theta_ijkl << std::endl;
	            }
	        }
	    }
	}
	
	std::cout << "------Torsion Angles------" << std::endl;
	for(int j = 0; j < s; j++){
	    for(int k = j+1; k < s; k++){
	        for(int i = 0; i < s; i++){
	            if(i == j || i == k)continue;
	            for(int l = i+1; l < s; l++){
	                if(l == j || l == k)continue;
	                double cos_tau_ijkl = tAngle(atoms[i].v, atoms[j].v, atoms[k].v, atoms[l].v);
	                std::cout << "cosτ(" << i+1 << j+1 << k+1 << l+1
	                      << ")= " << cos_tau_ijkl << std::endl;
	            }
	        }
	    }
	}
	
	std::cout << "------Center of Mass------" << std::endl;
	vec cmass = cMass(atoms);
	std::cout << "(" << cmass.x << ", " << cmass.y << ", " << cmass.z << ")" << std::endl;
	
	std::cout << "------Principle Moments of Inertia------" << std::endl;
	vec I = Minertia(atoms);
	std::cout << "Ia = " << I.x << "amu bohr2" << std::endl
	          << "Ib = " << I.y << "amu bohr2" << std::endl
	          << "Ic = " << I.z << "amu bohr2" << std::endl;
	          
	std::cout << "Ia = " << I.x * 0.280028 << "amu A2" << std::endl
	          << "Ib = " << I.y * 0.280028 << "amu A2" << std::endl
	          << "Ic = " << I.z * 0.280028 << "amu A2" << std::endl;

	std::cout << "Ia = " << I.x * 4.6498 * 1e-41 << "g cm2" << std::endl
	          << "Ib = " << I.y * 4.6498 * 1e-41 << "g cm2" << std::endl
	          << "Ic = " << I.z * 4.6498 * 1e-41 << "g cm2" << std::endl;
	          
	std::cout << "------Rotational Constants------" << std::endl;
    std::cout << "A = " << 60.19980/I.x << "cm^(-1)" << std::endl
	          << "B = " << 60.19980/I.y << "cm^(-1)" << std::endl
	          << "C = " << 60.19980/I.z << "cm^(-1)" << std::endl;
    std::cout << "A = " << 1804744.56/I.x << "MHz" << std::endl
	          << "B = " << 1804744.56/I.y << "MHz" << std::endl
	          << "C = " << 1804744.56/I.z << "MHz" << std::endl;

}