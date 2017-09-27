#include "non_interact.h"
#include <cmath>
#include <time.h>


Non_interact::Non_interact(double rho_max,vec n_list) {
    this->rho_max = rho_max;
    this->n_list = n_list;
}

void Non_interact::write_to_file(){
    ofstream outfile;
    mat energy = ones<vec>(3);
    mat known_values = vec({3.0, 7.0, 11.0});
    double eps = pow(10,-5);
    int iterations = 0;
    //int n = 400;
    double time;
    //int limit = n+2;

    outfile.open("../result.txt");
    outfile <<"rho_max = "<<rho_max <<endl;
    outfile << "Mesh points" << "    " << "iterations" << "    " << "time" << "            eigvec1  "<<"     eigvec2  "<<"       eigvec3 "<<endl;
    int i = 0;
    double n;
    //while(abs(energy(0) - known_values(0)) > eps &&
    //    abs(energy(1) - known_values(1)) > eps &&
    //    abs(energy(2) - known_values(2)) >eps) && (i<size(n_list)){
    while(i<size(n_list)[0]){
        n = n_list(i);
        Solve_SE_twoparticle(n, energy, iterations, time, rho_max);

        outfile << n << "               " << iterations << "           " << time << "          "<< energy(0) <<"       "<< energy(1) <<"       "<< energy(2) <<endl;
        cout << "n = "<< n<<endl;
        cout << "The eigenvalues are: " << energy(0) <<", "<< energy(1) <<" and "<< energy(2) << endl;
        //n +=1000;
        i +=1;

    }
   /* if (n==limit-1){
        cout << "Did not get four decimal points accuracy" << endl;
    }
    else{
        cout << n << " mesh points were necassary to get four decimal points." << endl;
    }
    */
    //cout << "rho_max: " << rho_max << "   " << "n: " << n << endl;

    outfile.close();

    }

void Non_interact::Solve_SE_twoparticle(int n, mat& energy, int& iterations, double& time, double rho_max){
    mat R, A;
    make_A(n, rho_max, A);

    //A.print("A");
    //cout << "results:"<<endl;
    //A.print("A");
    //R.print("R");

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
    time = (double) (finish_n - start_n)/(CLOCKS_PER_SEC );

    // Sorting eiegenvalues and eigenvectors
            mat eigenvalues = ones<vec>(n);
            mat eigenvectors = ones<mat>(n,n);

        for(int i=0;i<n;i++){
            eigenvalues(i) = A(i,i);
        }
  //      eigenvalues.print("lambda");

        uvec eigenindex =  sort_index(eigenvalues);

        eigenvalues = sort(eigenvalues);
        for(int j=0;j<n;j++){
            eigenvectors.col(j) = R.col(eigenindex(j));
        }

        //if(abs(min_element(eigenvalues) - energy(0))<= pow(10, -4)) cout << "its ok" << endl;

        energy(0) = eigenvalues(0);
        energy(1) = eigenvalues(1);
        energy(2) = eigenvalues(2);

}

void Non_interact::Jacobi(mat& A, mat& R, int n, int& iterations){
    double tol = 1e-10;
    int k=10; int l=10;
    double max_ = 10;
    int max_iterations = pow(10,5);
    iterations = 0;

    R = eye<mat>(n,n);
    /*for (int i =0;i<n;i++){
        for (int j =0;j<n;j++){
            if(i==j) R(i,j) = 1;
        }
    }
*/

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
    //mat rho = ones<vec>(n); // dimentionless radius
    double h = rho_max/(double) n;

    double rho;
    for(int i = 0; i<n;i++){
        rho = (i+1)*h;
        //rho(i) = rho;
        V(i) = rho*rho;
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

}

int Non_interact::test_eigensolver(){
    mat A = {{3,1},{1,3}};
    mat R;
    int iterations = 0;
//    A.print("Test matrix");
    Jacobi(A, R, 2, iterations);
//    A.print("Eigenvalues:");
//    R.print("Eigenvectors:");
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


