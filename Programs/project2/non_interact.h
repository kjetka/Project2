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

    Non_interact(vec n_list, int N_omega,  vec omega_list);


    double norm_off_diag(mat& A, int& k, int& l, int n);
    void write_to_file();
    void Jacobi_rot(mat& A, mat& R, int k, int l, int n);
    void Solve_SE_twoparticle(int n, mat& energy, int &iterations, double& time, double rho_max, double omega);
    int test_eigensolver();
    void make_A(int n, double max_rho, mat& A, double omega);
    void Jacobi(mat& A, mat& R, int n, int& iterations);
    int test_off_diagonal();

    void testing(double t);

    void armadillo_ref(int n, double rho_max, double &time_arma, double omega);
private:
    //double rho_max = 1;
    vec n_list = {1};
    vec omega_list = {1};
    int N_omega = 400;







};

#endif // NON_INTERACT_H



