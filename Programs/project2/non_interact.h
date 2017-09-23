#ifndef NON_INTERACT_H
#define NON_INTERACT_H
#include <armadillo>
#include <iostream>
#include <fstream>

using namespace arma;
using namespace std;

class Non_interact
{
public:

    Non_interact(double rho_max);


    void matrise(double rho_max, int n, mat& R, mat& A);
    double norm_off_diag(mat& A, int& k, int& l, int n);
    void write_to_file();
    mat Jacobi_rot(mat& A, mat& R, int k, int l, int n);
    void Jacobi_func(int N, mat& energy, int &iterations, double& time, double rho_max);
    void test_eigensolver();
private:
    int rho_max = 1;








};

#endif // NON_INTERACT_H



