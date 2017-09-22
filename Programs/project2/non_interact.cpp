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



double Non_interact::norm_off_diag(mat& A, int& k, int& l, int n){
    double max_a_kl =0;
    // int n = size(A)[0];

    for(int i = 0; i<n;i++){
        for(int j=i+1; j<n; j++){
           double aij = fabs(A(i,j));
            if(aij > max_a_kl){
                max_a_kl = aij;
                k = i; l = j;
            }
        }
    }
    return pow(max_a_kl,2);

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

mat Non_interact::Jacobi_rot(mat& A, mat& R, int k, int l, int n){
    //int n = size(A)[0];
    double s,c, t, tau, a_ll, a_kk, a_kl;
    double a_ik, a_il, r_ik, r_il;

     a_ll = A(l,l);  a_kk = A(k,k);
     a_kl = A(k,l);

     if(a_kl != 0.0){
        tau = (a_ll-a_kk)/(2.0*a_kl);
        if(tau>= 0){
             t = tau + sqrt(1+pow(tau,2));
        }
        else{
             t = -tau - sqrt(1+pow(tau,2));
        }
         c = 1.0/(    sqrt(1  +   pow(t,2))   );
         s = t*c;
     }
     else{
         c = 1.0;
         s = 0.0;
         }

     double cs = c*s; double cc = pow(c,2); double ss = pow(s,2);

     A(k,k) = cc * a_kk -2.0*cs*a_kl + ss*a_ll;
     A(l,l) = ss*a_kk + 2.0*cs*a_kl + cc*a_ll;
     A(k,l) = 0.0; A(l,l) = 0.0;
     cout <<A(k,k)<<endl;

     for(int i=0;i<n; i++){
         if (i!= k && i!= l){
             a_ik = A(i,k);
             a_il = A(i,l);
             A(i,k) = c*a_ik - s*a_il;
             A(k,i) = A(i,k);
             A(i,l) = c*a_il + s*a_ik;
             A(l,i) = A(i,l);
         }
         r_ik = R(i,k);
         r_il = R(i,l);
         R(i,k) = c*r_ik - s*r_il;
         R(i,l) = c*r_il + s*r_ik;
     }


            return A;
}


