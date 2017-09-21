#include <iostream>
#include "non_interact.h"
#include <armadillo>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <time.h>
#include <iomanip>

using namespace std;
using namespace arma;

Non_interact ni;


int main()
{

int N = 2;
mat V = ones<vec>(N); // potential
mat R = ones<mat>(N,N);
mat A = ni.matrise(V,N);

double tol = 1e-10;
int k; int l;
int iterations = 0;
double max_ = 10;

while(max_ > tol && iterations < 2){
    max_ = ni.norm_off_diag( A, k, l);
    A,R = ni.Jacobi_rot(A,R, k,l);
    iterations +=1;
    cout << iterations << endl;

}
A.print();


cout << "Problem: A,R will not return properly from Jacobi_rot!!!!!!!!!!!"<<endl;




cout << "A must be symmetric, should we test it?"<<endl;
cout << "minor adjust: remove finding n in every function..."<<endl;

cout << "Why t+ for tau >0 and t- for tau < 0?"<< endl;










}


