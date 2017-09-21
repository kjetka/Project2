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

double Non_interact::norm_off_diag(mat A, int k, int l){
    double max_akl =0;
    int n = size(A)[0];
    cout << "A must be symmetric, should we test it?"<<endl;

    for(int i = 0; i<n;i++){
        for(int j=i+1; j<n; j++){
           double aij = fabs(A(i,j));
            if(aij > max_akl){
                max_akl = aij;
                k = i; l = j;
            }
        }
    }
    return max_akl;
    /*
    double off_A_verdi=0;
    for(int i=0;i<n;i++){
        for(int j=0; j<n;j++){
            if(j!=i)
                off_A_verdi += A(i,j);


        }
    }
    off_A_verdi = sqrt(off_A_verdi);

    return off_A_verdi;
    */

}
