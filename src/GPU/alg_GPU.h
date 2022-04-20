#include<iostream>

/**
 forward declarations of conjugate gradient algos for double and float fp numbers.
*/

/** conjugate gradient for double */
double cg(int *I, int *J, double *val, double *x, double *rhs, const int N, const double tol, const int max_iter, int &nb_iter);

/** conjugate gradient for float */
float cg(int *I, int *J, float *val, float *x, float *rhs, const int N, const float tol, const int max_iter, int &nb_iter);


/** print the cuda archtecture and the number of multi-processors of the GPU */
void infos(void);

