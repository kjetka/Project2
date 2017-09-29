#include "non_interact.h"
#include <cmath>
#include <time.h>
#include <armadillo>

using namespace arma;
using namespace std;

Non_interact::Non_interact(double rho_max, vec n_list, double omega){
    this->rho_max = rho_max;
    this->n_list = n_list;
    this->omega = omega;
}

void Non_interact::write_to_file(){
    mat energy = ones<vec>(3);
    mat known_values = vec({3.0, 7.0, 11.0});
    int iterations = 0;

    double time_jacobi;
    string str = to_string(omega);
    str.resize(4);
    string filename = string("../outfiles/")+ string("result_omega_") + str + string(".txt");

    ofstream outfile;


    if(omega == 0) {
        ofstream outfile1;
        ofstream outfile2;


        double time_arma;
        string  filename1 = "../outfiles/result_eigval.tex";
        string filename2 = "../outfiles/result_time.tex";

        outfile1.open(filename1);
        outfile2.open(filename2);
        outfile1 << "N   &   $\\lambda_1$  &  $\\lambda_2$  &  $\\lambda_3$    "<<endl;
        outfile1 << "\\hline"<< endl;
        outfile2 <<"N       &   # Similarity transforms   &  time Jacobi (s)      & time Armadillo (ms)            "<<endl;
        outfile2 << "\\hline"<< endl;

        for (int i=0; i<size(n_list)[0]; i++){
            iterations = 0;
            time_arma = 0;
            time_jacobi =0;
            double n = n_list(i);
            Solve_SE_twoparticle(n, energy, iterations, time_jacobi, rho_max);
            armadillo_ref(n, rho_max, time_arma);


            outfile1 << n << "     &          " <<  energy(0) <<"   &    "<< energy(1) <<"   &    "<< energy(2) <<endl;
            outfile2 << n<< "      &   "   << iterations  <<  "      &   " << time_jacobi << "      &   " << time_arma*1e3 << endl;

            cout << "time arma:  "<<time_arma <<endl;
            cout << "time jacobi:  "<<time_jacobi <<endl;

        }

        outfile1.close();
        outfile2.close();




    }


    //-------------------------------------

    //-------------------------------------

else{




    outfile.open(filename);
    outfile <<"rho_max = "<<rho_max <<endl;

    }
    outfile << "N   &      eigvec1    &       eigvec2    &      eigvec3 "<<endl;
    outfile << "\\hline"<<endl;

    for (int i=0; i<size(n_list)[0]; i++){

        double n = n_list(i);
        Solve_SE_twoparticle(n, energy, iterations, time_jacobi, rho_max);


        outfile << n << "               "    << "          "<< energy(0) <<"       "<< energy(1) <<"       "<< energy(2) <<endl;
    }

    outfile.close();
    }





void Non_interact::Solve_SE_twoparticle(int n, mat& energy, int& iterations, double& time_jacobi, double rho_max){
    mat R, A;
    make_A(n, rho_max, A);

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

  //      eigenvalues.print("lambda");

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
    int statement = 1;//test_off_diagonal();
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

void Non_interact::make_A(int n, double rho_max, mat& A){

    A = zeros<mat>(n,n);
    mat V = ones<vec>(n); // potential
    double h = rho_max/(double) n;
    double rho;

    for(int i = 0; i<n;i++){
            rho = (i+1)*h;
        if(omega == 0){
            //rho(i) = rho;
            V(i) = rho*rho;
        }
        else{
            V(i) = omega*omega*rho*rho + 1/rho; //OBS! Blir rho veldig liten? evt veldig stor?
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
    return max_a_kl;

}

void Non_interact::Jacobi_rot(mat& A, mat& R, int k, int l, int n){
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

    for(int i = 0; i<n;i++){
        for(int j=i+1; j<n; j++){
           double aij = fabs(A(i,j));
            if(aij > max_a_kl){
                max_a_kl = aij;
            }
        }
    }
    if (abs(max_a_kl- 3.0) <= pow(10,-4)){
        return 1;
    }
    else{
        return 0;
    }
}

void Non_interact::armadillo_ref(int n, double rho_max, double& time_arma){

       mat A = zeros<mat>(n,n);
       make_A(n, rho_max, A);
       clock_t start_arma, finish_arma;
       start_arma = clock();

       vec eigval = eig_sym(A);

       finish_arma = clock();
       time_arma = (double) (finish_arma - start_arma)/(CLOCKS_PER_SEC );
}




