#include <iostream>
#include "non_interact.h"
#include <armadillo>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <time.h>
#include <iomanip>
#include <iostream>
#include <fstream>

using namespace std;
using namespace arma;

int main(){

    Non_interact* ni;


    double rho_max;

    for(int rhos=1;rhos<20;rhos++){
        rho_max = rhos*1.0;
        ni = new Non_interact(rho_max);
        ni->write_to_file();
    }

//    ni->test_eigensolver();
    //ni->Jacobi_func(3);


cout << "A must be symmetric, should we test it?"<<endl;









}


