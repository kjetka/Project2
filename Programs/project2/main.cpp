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
    clock_t start1, finish1;
    start1 = clock();

   // for(int rhos=1;rhos<20;rhos++){
        rho_max = 50.0; //rhos*1.0;
        ni = new Non_interact(rho_max);
        ni->write_to_file();


    //}

        finish1 = clock();
        double tid = (double) (finish1 - start1)/(CLOCKS_PER_SEC );
cout << "Time entire crap "<<tid<<endl;
cout << "A must be symmetric, should we test it?"<<endl;









}


