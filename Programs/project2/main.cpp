#include <iostream>
#include "Non_interact.h"
#include <armadillo>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <time.h>
#include <iomanip>


using namespace std;
using namespace arma;

int main()
{
Non_interact ni;

cout << ni.test( 2.0)<<endl;
mat V = ones<vec>(3);
mat A = ni.matrise(V,3);
A.print();
















}


