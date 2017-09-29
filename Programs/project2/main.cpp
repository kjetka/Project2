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
    Non_interact* nipluss;

    vec n_list= {10, 100, 200,400, 600};
    vec omega_list =    {0,     0.01,   0.5,    1,      5};
    vec rho_max_list =  {4.79,  10,     10,     4.79,   4.79};

    if (size(rho_max_list)[0] != size(omega_list)[0]) {
        exit(10);
    }

    int stopp = size(omega_list)[0];
    for(int i = 0;i<stopp;i=i+2){
        ni = new Non_interact(rho_max_list(i), n_list, omega_list(i)); //for the non-interacting case omega = 0
        ni->write_to_file();
        if (i<size(omega_list)[0]){
        nipluss = new Non_interact(rho_max_list(i+1), n_list, omega_list(i+1));
        nipluss->write_to_file();
        }

        cout<< "Rho_max = "<< fixed<< rho_max_list(i) << endl;
     cout << "---------------------" <<endl;
        cout << "rho_max = "<< rho_max_list(i)<<endl;
        for (int j=0; j < stopp; j++){
            double n = n_list(j);
            cout    << "n = "<< n<<endl;
           cout<< "h = "<<rho_max_list(i)/n << endl;
    }
}










}


