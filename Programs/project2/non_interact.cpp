#include "non_interact.h"
#include <cmath>


Non_interact::Non_interact(double rho_max) {
    this->rho_max = rho_max;
}

void Non_interact::write_to_file(){
    ofstream outfile;
    mat energy = ones<vec>(3);
    mat known_values = vec({3.0, 7.0, 11.0});
    double eps = pow(10,-5);
    int iterations = 0;
    int N = 3;
    double time;

    outfile.open("result.txt");
    outfile << "Mesh points" << "    " << "iterations" << "    " << "time" << endl;

    while(abs(energy(0) - known_values(0)) > eps &&
        abs(energy(1) - known_values(1)) > eps &&
        abs(energy(2) - known_values(2)) >eps && N < 20){
        Jacobi_func(N, energy, iterations, time, rho_max);
        outfile << N << "    " << iterations << "    " << time << endl;
        N +=1;
}
N = N-1;
cout << N << " mesh points were necassary to get four decimal points." << endl;
cout << "The eigenvalues are: " << energy(0) <<", "<< energy(1) <<" and "<< energy(2) << endl;
outfile.close();

}

void Non_interact::Jacobi_func(int N, mat& energy, int& iterations, double& time, double rho_max){
    double tol = 1e-6;
    int k=10; int l=10;
    double max_ = 10;
    int max_iterations = N*100;

    mat R, A;
    matrise(rho_max, N, R, A);

    while(max_ > tol && iterations < max_iterations){
        max_ = norm_off_diag( A, k, l, N);
        A = Jacobi_rot(A,R, k,l, N);
        iterations +=1;
    }
    /*
    mat printA = A;
    for(int i = 0; i<N;i++){
        for(int j = 0; j<N;j++){
            if(i != j && printA(i,j) <= tol){
                printA(i,j)= 0;
           }
        }
    }

   if(N==50){
       printA.print("A");
       cout << "iterations: " << iterations << endl;
   }
*/
    iterations = iterations-1;

    //cout << "results:"<<endl;
    //A.print("A");
    //R.print("R");

    mat eigenvalues = ones<vec>(N);
    mat eigenvectors = ones<mat>(N,N);

        for(int i=0;i<N;i++){
            eigenvalues(i) = A(i,i);
        }
        eigenvalues.print("lambda");

        uvec eigenindex =  sort_index(eigenvalues);

        eigenvalues = sort(eigenvalues);
        for(int j=0;j<N;j++){
            eigenvectors.col(j) = R.col(eigenindex(j));
        }
        energy(0) = eigenvalues(0);
        energy(1) = eigenvalues(1);
        energy(2) = eigenvalues(2);
        time = 8.0;
}

void Non_interact::matrise(double rho_max, int n, mat& R, mat &A){
    A = zeros<mat>(n,n);

    mat V = ones<vec>(n); // potential
    mat rho = ones<vec>(n); // dimentionless radius
    double h = rho_max/(double) n;

    for(int i = 0; i<n;i++){
        rho(i) = i*h;
        V(i) = rho(i)*rho(i);
    }


    for (int i=0; i<n; i++) {
        A(i,i)=2.0/(double)(h*h) + V[i];
        if (i!=n-1) A(i,i+1) = -1.0/(double)(h*h);
        if (i!=n-1) A(i+1,i) = -1.0/(double)(h*h);
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
     //cout <<A(k,k)<<endl;

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

void Non_interact::test_eigensolver(){
    mat B = {{3,1},{1,3}};
    int iterations;
    double time;
    double rho_max;
    mat eigenvalues;
    B.print();
    Jacobi_func(2, eigenvalues, iterations , time, rho_max);
    eigenvalues.print("Test eigenvalues");
}

