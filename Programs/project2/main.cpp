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

    vec n_list= {10,100, 200,400};
    int N_omega = 400;
    vec omega_list =    {0,     0.01,   0.5,    1,      5};
    //vec rho_max_list =  {4.79,  10,     10,     10,     10};


        ni = new Non_interact(n_list,N_omega,  omega_list); //for the non-interacting case omega = 0
        ni->write_to_file();






}


