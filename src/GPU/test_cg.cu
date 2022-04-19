#include "alg_GPU.h"
//#include "cg.h"


#include<iomanip>

#include <stdlib.h>
#include <time.h>

#include<type_traits>

#include "../alg_utils.h"


/* genTridiag: generate a random tridiagonal symmetric matrix , zero based indices */
template <typename T>
void genTridiag(int *I, int *J, T *val, int N, int nz)
{
srand(time(NULL));

  I[0] = 0, J[0] = 0, J[1] = 1;
  val[0] = (T) rand() / RAND_MAX + 10.0f;
  val[1] = (T) rand() / RAND_MAX;
  int start;

  for (int i = 1; i < N; i++) 
  {
    if (i > 1) { I[i] = I[i - 1] + 3; } else { I[1] = 2; }

    start = (i - 1) * 3 + 2;
    J[start] = i - 1;
    J[start + 1] = i;

    if (i < N - 1) { J[start + 2] = i + 1; }

    val[start] = val[start - 1];
    val[start + 1] = (T) rand() / RAND_MAX + 10.0f;

    if (i < N - 1) { val[start + 2] = (T) rand() / RAND_MAX; }
  }

  I[N] = nz;
}


template <typename T>
T check_sol(int *I,int *J, T *val, T *x,int N, T *rhs)
{
T result(0);
T *y;

y = new T[N];

alg::multCSR_MatVect<T>(I,J,val,x,N,y);

for(int k =0;k<N;k++) 
	{ result += alg::sq<T>( y[k] - rhs[k]); }
delete [] y;

return sqrt(result);
}

int main(void)
{
infos();

const int N =100000;
int nz = (N-2)*3 + 4;
int *I, *J ;
I = new int[N+1];
J = new int[nz];

double tol,res;

decltype(tol) *val, *x, *rhs;

val = new decltype(tol)[nz];
x = new decltype(tol)[N];
rhs = new decltype(tol)[N];

genTridiag<decltype(tol)>(I, J, val, N, nz);

for (int i = 0; i<N; i++)
	{
	rhs[i] = 1.0;
	x[i] = rhs[i]; // initial guess
	}

tol = 1e-6;
int max_iter(100);
int nb_iter(0);

res = cg(I,J,val,x,rhs,N,tol,max_iter,nb_iter);

//res = cg<decltype(tol)>(I,J,val,x,rhs,N,tol,max_iter,nb_iter);
std::cout << "nb iter = " << nb_iter << "; residu = " << res << std::endl;

std::cout << "check solution returns : " << check_sol(I,J,val,x,N,rhs) <<std::endl;

delete[] x;

delete[] I;
delete[] J;
delete[] val;
delete[] rhs;

std::cout <<"CUDA error: " << cudaGetErrorString(cudaGetLastError()) << std::endl;
return 0;
}
