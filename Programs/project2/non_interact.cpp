#include "non_interact.h"
#include <cmath>
#include <armadillo>



Non_interact::Non_interact() {
}


mat Non_interact::matrise(mat V, int n)
{

    mat A_ = zeros<mat>(n,n);
    for (int i=0; i<n; i++) {
        A_(i,i)=2 + V[i];
        if (i!=n-1) A_(i,i+1) = 1;
        if (i!=n-1) A_(i+1,i) = 1;
}
    return A_;

}

double Non_interact::off(mat A){
    double verdi=0;
    int n = size(A)[0];
    for(int i=0;i<n;i++){
        for(int j=0; j<n;j++){
            if(j!=i)
                verdi += A(i,j);

        }
    }

    return verdi;

}
