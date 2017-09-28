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
    vec n_list= {10, 100,200,400};
    vec omega_list = {0, 0.01, 0.5, 1.0, 5.0};

   // for(int rhos=1;rhos<20;rhos++){

        cout<< "Rho_max = "<< fixed<< rho_max << endl;

        for(int i = 0;i<5;i++){
            ni = new Non_interact(rho_max, n_list, omega_list(i)); //for the non-interacting case omega = 0
            ni->write_to_file();
        }



    //}
/*
        finish1 = clock();
        double tid = (double) (finish1 - start1)/(CLOCKS_PER_SEC );


cout << "Time entire crap "<<tid<<endl;










}


