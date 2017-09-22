#include "non_interact.h"
#include <cmath>
#include <armadillo>
#include <iostream>


Non_interact::Non_interact() {
    //this->N = N;
}

void Non_interact::Jacobi_func(int N){
    double tol = 1e-10;
    int k=10; int l=10;
    int iterations = 0;
    double max_ = 10;
    int max_interations = N*10;
    mat V = ones<vec>(N); // potential

    mat R, A;
    matrise(V, N, R, A);

    while(max_ > tol && iterations < max_interations){
        max_ = norm_off_diag( A, k, l, N);
        A = Jacobi_rot(A,R, k,l, N);
        iterations +=1;
    }

    cout << "Iterations: "<<iterations << endl;

    cout << "results:"<<endl;
    A.print("A");
    R.print("R");

    mat eigenvalues = ones<vec>(N);
    mat eigenvectors = ones<mat>(N,N);
    mat eigenindex = ones<vec>(N);

        for(int i=0;i<N;i++){
            eigenvalues(i) = A(i,i);
        }
        eigenindex = sort_index(eigenvalues);
        eigenvalues = sort(eigenvalues);
        for(int j=0;j<N;j++){
            eigenvectors.col(j) = R.col(eigenindex(j));
        }
}
/*
void Non_interact::print_to_file(int iterations, int N, double time){
}
*/
void Non_interact::matrise(mat V, int n, mat& R, mat &A){
    A = zeros<mat>(n,n);
    for (int i=0; i<n; i++) {
        A(i,i)=2 + V[i];
        if (i!=n-1) A(i,i+1) = -1;
        if (i!=n-1) A(i+1,i) = -1;
}
    R = zeros<mat>(n,n);
    for (int i =0;i<n;i++){
        for (int j =0;j<n;j++){
            if(i==j) R(i,j) = 1;
        }
    }

    return;
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
             t = 1.0/(tau + sqrt(1+tau*tau));
        }
        else{
            t = -1.0/(- tau + sqrt(1+tau*tau));
        }

         c = 1.0/(    sqrt(1  +  t*t)   );
         s = t*c;
     }
     else{
         c = 1.0;
         s = 0.0;
         }

     double cs = c*s; double cc = pow(c,2); double ss = pow(s,2);

     A(k,k) = cc * a_kk -2.0*cs*a_kl + ss*a_ll;
     A(l,l) = ss*a_kk + 2.0*cs*a_kl + cc*a_ll;
     A(k,l) = 0.0; A(l,k) = 0.0;
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


