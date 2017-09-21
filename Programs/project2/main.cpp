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
mat V = ones<vec>(N);
mat A = ni.matrise(V,N);
double tall = ni.off(A);
cout <<tall<<endl;
A.print();













}


