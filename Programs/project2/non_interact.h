#ifndef NON_INTERACT_H
#define NON_INTERACT_H
#include <armadillo>
using namespace arma;

class Non_interact
{
public:
    Non_interact();
    double test(double masse);

    mat matrise(mat V, int n);
};

#endif // NON_INTERACT_H
