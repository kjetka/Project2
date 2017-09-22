#ifndef NON_INTERACT_H
#define NON_INTERACT_H
#include <armadillo>
using namespace arma;

class Non_interact
{
public:
    Non_interact();
    void matrise(mat V, int n, mat& R, mat& A);
    double norm_off_diag(mat& A, int& k, int& l, int n);

    mat Jacobi_rot(mat& A, mat& R, int k, int l, int n);





private:









};

#endif // NON_INTERACT_H



