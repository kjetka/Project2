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


    clock_t start1, finish1;
    start1 = clock();
    vec n_list= {10, 100, 200,400};
    vec omega_list = {0,0.01, 0.5,1,5};
    vec rho_max_list = {4.79,10, 10,4.79,4.79};

    if (size(rho_max_list)[0] != size(omega_list)[0]) {
        exit(10);
    }

   // for(int rhos=1;rhos<20;rhos++){

        int stopp = size(omega_list)[0];
        for(int i = 0;i<stopp;i++){
            ni = new Non_interact(rho_max_list(i), n_list, omega_list(i)); //for the non-interacting case omega = 0
            ni->write_to_file();
            cout<< "Rho_max = "<< fixed<< rho_max_list(i) << endl;

        }



    //}

        finish1 = clock();
        double tid = (double) (finish1 - start1)/(CLOCKS_PER_SEC );


cout << "Time entire crap "<<tid<<endl;










}


