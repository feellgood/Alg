#include<iomanip>
#include <stdlib.h>
#include <time.h>
#include<type_traits>

#include "alg_GPU.h" // headers for the declaration of the functions of the static library alg_GPU
#include "../alg_utils.h"


/* genTridiag: generate a random tridiagonal symmetric matrix , zero based indices */
template <typename T>
void genTridiag(alg::CSR_mat<T> &A)
{
srand(time(NULL));

A.I[0] = 0, A.J[0] = 0, A.J[1] = 1;
A.val[0] = (T) rand() / RAND_MAX + 10.0f;
A.val[1] = (T) rand() / RAND_MAX;
int start;

for (int i = 1; i < A.N; i++) 
  {
    if (i > 1) { A.I[i] = A.I[i - 1] + 3; } else { A.I[1] = 2; }

    start = (i - 1) * 3 + 2;
    A.J[start] = i - 1;
    A.J[start + 1] = i;

    if (i < A.N - 1) { A.J[start + 2] = i + 1; }

    A.val[start] = A.val[start - 1];
    A.val[start + 1] = (T) rand() / RAND_MAX + 10.0f;

    if (i < A.N - 1) { A.val[start + 2] = (T) rand() / RAND_MAX; }
  }
}


template <typename T>
T check_sol(alg::CSR_mat<T> const& A,T *x,T *rhs)
{
T result(0);

//T *y;
const int DIM_y = A.N;

T y[DIM_y];
//y = new T[DIM_y];


for(int k =0;k<DIM_y;k++) {y[k]= 0;} //initialization to zero

alg::multCSR_MatVect<T>(A,x,y);

for(int k =0;k<DIM_y;k++) 
	{ result += alg::sq<T>( y[k] - rhs[k]); }
//delete [] y;

return sqrt(result);
}

int main(void)
{
infos();

double tol,res;

const int N =100000;

const int nb_coeff = (N-2)*3 + 4;// nb_coeff of a tri_diag sparse mat

alg::CSR_mat<decltype(tol)> A(N,nb_coeff);

genTridiag<decltype(tol)>(A);
std::cout << "random tri-diagonal sparse matrix filled.\n";

decltype(tol) *x, *rhs;
x = new decltype(tol)[N];
rhs = new decltype(tol)[N];

for (int i = 0; i<N; i++)
	{
	rhs[i] = 1.0;
	x[i] = rhs[i]; // initial guess
	}

tol = 1e-6;
int max_iter(100);
int nb_iter(0);

std::cout << "now starting gradient conjugate on GPU... ";
res = cg(A,x,rhs,tol,max_iter,nb_iter);
std::cout << "\t ... job done on GPU.\n";

std::cout << "nb iter = " << nb_iter << "; residu = " << res << "\tcheck solution returns : " << check_sol(A,x,rhs) <<std::endl;

delete[] x;
delete[] rhs;

std::cout <<"CUDA error: " << cudaGetErrorString(cudaGetLastError()) << std::endl;
return 0;
}
