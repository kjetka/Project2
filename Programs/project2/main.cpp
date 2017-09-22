#include <iostream>
#include "non_interact.h"
#include <armadillo>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <time.h>
#include <iomanip>
#include <iostream>

using namespace std;
using namespace arma;


ofstream outfile;

int main(){

    Non_interact* ni;

    //ni = new Non_interact(2);


    mat B = ones<vec>(3);
    B(1) = 3; B(2) = 0;
    B.print();
    cout<< "-------"<<endl;
    cout << sort_index(B)<<endl;
    /*
    for(int i=0;i<3;i++){
        int j = sort_index(B)[i];
        cout<< B(j)<<endl;
    }*/

//ni->Jacobi_func(2);


cout << "A must be symmetric, should we test it?"<<endl;
cout << "Why t+ for tau >0 and t- for tau < 0?"<< endl;









}


