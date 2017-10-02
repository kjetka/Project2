#include "non_interact.h"
#include <cmath>
#include <time.h>
#include <armadillo>
#include <fstream>
#include <iostream>
#include "non_interact.h"
#include <cstdlib>
#include <iomanip>
#include <iostream>

using namespace arma;
using namespace std;

Non_interact::Non_interact(vec n_list,int N_omega, vec omega_list){
    this->n_list = n_list;
    this-> N_omega = N_omega;
    this->omega_list = omega_list;
}


void Non_interact::write_to_file(){
    mat energy = ones<vec>(3);
    mat known_values = vec({3.0, 7.0, 11.0});
    int iterations = 0;

    /*
    double bee = 1.44;
    double m = 0.51099*pow(10,6);
    double hbarc2 =pow(1240.0/(2*M_PI),2);
    double alfa =   hbarc2/(m*bee);
    double E = hbarc2/(m*alfa*alfa);
    */

    double time_jacobi;
    string filename = string("../../outfiles/results_omega.txt");



    ofstream outfile;

    outfile.open(filename);

    double omega;
    double rho_max ;
    int stopp = size(omega_list)[0];
    outfile << "   $\\omega_r$  &  $\\lambda_0 $   &       $\\lambda_1$  &      $\\lambda_2$ \\\\"<<endl;
    outfile << "\\hline\\"<<endl;

    for(int i = 0;i<stopp;i++){
    omega = omega_list(i);
    if(omega == 0) {
        ofstream outfile1;
        ofstream outfile2;
        rho_max = 4.79;

        double time_arma;
        string  filename1 = "../../outfiles/result_eigval.txt";
        string filename2 = "../../outfiles/result_time.txt";

        outfile1.open(filename1);
        outfile2.open(filename2);
        outfile1 << "N   &   $\\lambda_1$  &  $\\lambda_2$  &  $\\lambda_3$  \\\\  "<<endl;
        outfile1 << "\\hline"<< endl;
        outfile2 <<"N       &   Transforms   &  time Jacobi (s)      & time Armadillo (ms)     \\\\       "<<endl;
        outfile2 << "\\hline"<< endl;
        outfile2 << setprecision(3);
        for (int i=0; i<size(n_list)[0]; i++){
            iterations = 0;
            time_arma = 0;
            time_jacobi =0;
            double n = n_list(i);
            Solve_SE_twoparticle(n, energy, iterations, time_jacobi, rho_max, omega);
            armadillo_ref(n, rho_max, time_arma, omega);


            outfile1 <<n << "     &          " <<  energy(0) <<"   &    "<< energy(1) <<"   &    "<< energy(2)<<"\\\\" <<endl;

            outfile2 <<defaultfloat<< n<< "      &   "   << iterations  <<  "      &   " << fixed<< time_jacobi << "      &   " << time_arma*1e3<<"\\\\" << endl;


            cout << "time arma:  "<<time_arma <<endl;
            cout << "time jacobi:  "<<time_jacobi <<endl;

        }

        outfile1.close();
        outfile2.close();




    }


    else{

        if (omega<0.1) rho_max = 50;
        else rho_max = 7.9; //rho_max = 4.79;


        cout << "Rho = "<<rho_max <<endl;
        double n = N_omega;
        Solve_SE_twoparticle(n, energy, iterations, time_jacobi, rho_max, omega);


        outfile << defaultfloat<<omega << "       &        " << setprecision(3) <<fixed << energy(0) <<"    &   "<< energy(1) <<"   &    "<< energy(2)<<"\\\\" <<endl;


    }

}
    outfile.close();

}

void Non_interact::Solve_SE_twoparticle(int n, mat& energy, int& iterations, double& time_jacobi, double rho_max, double omega){
    mat R, A;
    make_A(n, rho_max, A, omega);

    // Test if Jacobi works
    int statement = test_eigensolver();

    clock_t start_n, finish_n;
    start_n = clock();

    if (statement == 1){
        Jacobi(A, R, n, iterations);
    }
        else{
        cout << "Error! Something is wrong in Jacobi. Getting wrong eigenvalues." << endl;
        exit(2);
    }
    finish_n = clock();
    time_jacobi = (double) (finish_n - start_n)/double((CLOCKS_PER_SEC ));

    // Sorting eiegenvalues and eigenvectors
            mat eigenvalues = ones<vec>(n);
            //mat eigenvectors = ones<mat>(n,n);

        for(int i=0;i<n;i++){
            eigenvalues(i) = A(i,i);
        }


        uvec eigenindex =  sort_index(eigenvalues);

        eigenvalues = sort(eigenvalues);
        /*
        for(int j=0;j<n;j++){
            eigenvectors.col(j) = R.col(eigenindex(j));
        }

        */
        energy(0) = eigenvalues(0);
        energy(1) = eigenvalues(1);
        energy(2) = eigenvalues(2);

}

void Non_interact::Jacobi(mat& A, mat& R, int n, int& iterations){
    double tol = 1e-10;
    int k=10; int l=10;
    double max_ = 10;
    int max_iterations = pow(10,6);
    iterations = 0;

    R = eye<mat>(n,n);
    // Test: Is it finding the maximum value?
    int statement = test_off_diagonal();
    if (statement == 1){
        while(max_ > tol && iterations < max_iterations){
            max_ = norm_off_diag( A, k, l, n);
            Jacobi_rot(A,R, k,l, n);
            iterations +=1;
        }
    }
        else{
        cout << "Error! Something is wrong in norm_off_diagonal. Can't find the maximum value." << endl;
        exit(3);
    }
}

void Non_interact::make_A(int n, double rho_max, mat& A, double omega){

    A = zeros<mat>(n,n);
    mat V = ones<vec>(n); // potential
    double h = rho_max/(double) n;
    double rho;

    for(int i = 0; i<n;i++){
            rho = (i+1)*h;
        if(omega == 0){
            V(i) = rho*rho;
        }
        else{
            V(i) = omega*omega*rho*rho + 1/rho;
        }
    }

    for (int i=0; i<n; i++) {
        double hh = h*h;
        A(i,i)= 2.0/(hh) + V[i];
        if (i!=n-1) {
            A(i,i+1) = -1.0/(hh);
            A(i+1,i) = -1.0/(hh);
        }
}
    return;
}

double Non_interact::norm_off_diag(mat& A, int& k, int& l, int n){
    double max_a_kl =0;

    for(int i = 0; i<n;i++){
        for(int j=i+1; j<n; j++){
           double aij = fabs(A(i,j));
            if(aij > max_a_kl){
                max_a_kl = aij;
                k = i; l = j;
            }
        }
    }
    return max_a_kl;

}

void Non_interact::Jacobi_rot(mat& A, mat& R, int k, int l, int n){
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

}

int Non_interact::test_eigensolver(){
    mat A = {{3,1},{1,3}};
    mat R;
    int iterations = 0;
    Jacobi(A, R, 2, iterations);
    if (abs(A(0,0)- 2.0) <= pow(10,-4) && abs(A(1,1)- 4.0) <= pow(10,-4)){
        return 1;
    }
    else{
        return 0;
    }
}

int Non_interact::test_off_diagonal(){
    mat A = {{1, 3, 1},{2, 1, 0.5},{1.5, 6, 2}};

    double max_a_kl =0;
    int n = size(A)[0];
    int k =0; int l = 0;
    norm_off_diag(A, k, l,  n);
    if (abs(  A(k,l)- 3.0) <= pow(10,-8)){
        return 1;
    }
    else{
        return 0;
    }
}

void Non_interact::armadillo_ref(int n, double rho_max, double& time_arma, double omega){

       mat A = zeros<mat>(n,n);
       make_A(n, rho_max, A, omega);
       clock_t start_arma, finish_arma;
       start_arma = clock();

       vec eigval = eig_sym(A);

       finish_arma = clock();
       time_arma = (double) (finish_arma - start_arma)/(CLOCKS_PER_SEC );
}




