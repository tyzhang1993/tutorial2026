#include"box.h"
#include<algorithm>

vec Eigenv(double A[3][3]){
    const double convergence = 1e-12;
    
    while(true){
        for(int i = 0; i < 3; i++){
            for(int j = i + 1; j < 3; j++){
                if(std::abs(A[i][j]) < 1e-15)continue; //  in case the denominator equals 0.
                double tau = (A[i][i] - A[j][j])/(2 * A[i][j]);
                double t;
                if(tau >= 0){
                    t = 1 / (tau + std::sqrt(1 + tau*tau)); 
                }
                else{
                    t = -1 / (-tau + std::sqrt(1 + tau*tau));
                }
                double c = 1 / std::sqrt(1 + t*t);
                double s = t * c;
                
                A[i][i] = A[i][i] - t*A[i][j];
                A[j][j] = A[j][j] + t*A[i][j];
                A[i][j] = 0.0;
                A[j][i] = 0.0;
                
                for(int k = 0; k < 3; k++){
                    if(k == i || k == j)continue;
                    double oldik = A[i][k];
                    double oldjk = A[j][k];
                    A[i][k] = c*oldik - s*oldjk;
                    A[j][k] = s*oldik + c*oldjk;
                    A[k][i] = A[i][k];
                    A[k][j] = A[j][k];
                }
            }
        }
        
        double max = 0.0;
        for(int i = 0; i < 3; i++){
            for(int j = i + 1; j < 3; j++){
                if(std::abs(A[i][j]) > max){
                    max = std::abs(A[i][j]);
                }
                else{
                    continue;
                }
            }
        }
        
        if(max < convergence)break;
    }
    
    double Ia = A[0][0];
    double Ib = A[1][1];
    double Ic = A[2][2];
    
    if(Ia > Ib)std::swap(Ia, Ib);
    if(Ia > Ic)std::swap(Ia, Ic);
    if(Ib > Ic)std::swap(Ib, Ic);
    
    return {
        Ia, Ib, Ic
    };
}