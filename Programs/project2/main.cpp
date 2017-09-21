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

int N = 4;
mat V = ones<vec>(N); // potential
mat R = zeros<mat>(N,N);
mat A = ni.matrise(V,N);
int k; int l;
double max_ = ni.norm_off_diag( A, k, l);
cout << "minor adjust: remove finding n in every function..."<<endl;











}


