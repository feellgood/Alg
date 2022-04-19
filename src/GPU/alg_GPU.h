#include<iostream>

/*
extern "C" 
{
*/

double cg(int *I, int *J, double *val, double *x, double *rhs, const int N, const double tol, const int max_iter, int &nb_iter);
//float cg(int *I, int *J, float *val, float *x, float *rhs, const int N, const float tol, const int max_iter, int &nb_iter);

//template <typename T> T cg(int *I, int *J, T *val, T *x, T *rhs, const int N, const T tol, const int max_iter, int &nb_iter);


void infos(void);

//}
