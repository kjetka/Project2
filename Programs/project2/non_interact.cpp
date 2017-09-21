#include "non_interact.h"
#include <cmath>
#include <iostream>
#include <armadillo>
#include <cstdlib>
#include <fstream>
#include <time.h>
#include <iomanip>


Non_interact::Non_interact(){

}

double Non_interact::test(double masse){
    return pow(masse,2);
}

mat Non_interact::matrise(mat V, int n)
{
    mat A_ = zeros<mat>(n,n);
    for (int i=0; i<n; i++) {
        A_(i,i)=2 + V[i];
        if (i!=n-1) A_(i,i+1) = 1;
        if (i!=n-1) A_(i+1,i) = 1;
}
    return A_;

}
