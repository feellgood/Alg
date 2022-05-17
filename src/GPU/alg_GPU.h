#include<iostream>
#include "../alg_utils.h"
/**
 forward declarations of conjugate gradient algos for double and float floating point numbers.
*/

/** conjugate gradient for double */
double cg(alg::CSR_mat<double> const& A, double *x, double *rhs, const double tol, const int max_iter, int &nb_iter);

/** conjugate gradient for float */
float cg(alg::CSR_mat<float> const& A, float *x, float *rhs, const float tol, const int max_iter, int &nb_iter);


/** print the cuda architecture and the number of multi-processors of the GPU */
void infos(void);


/** print last CUDA error */
void status(void);
