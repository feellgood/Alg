#include "cg.h"

#include<iostream>
#include<iomanip>

#include<type_traits>

/* genTridiag: generate a random tridiagonal symmetric matrix */
template <typename T>
void genTridiag(int *I, int *J, T *val, int N, int nz)
{
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

int main(void)
{
const int N =10;
int nz = (N-2)*3 + 4;
int *I, *J ;
I = new int[N+1];
J = new int[nz];

double tol;

decltype(tol) *val, *x, *rhs;

val = new decltype(tol)[nz];
x = new decltype(tol)[N];
rhs = new decltype(tol)[N];

genTridiag<decltype(tol)>(I, J, val, N, nz);
//for(int i=0;i<nz;i++)	{std::cout<< "val[" << i << "]: " << val[i] <<std::endl;}


for (int i = 0; i<N; i++)
	{
	x[i] = i/10.f;
	rhs[i] = 1.0;
	std::cout<< "x[" << i << "]: " << x[i] <<std::endl;
	}

tol = 1e-6;
int max_iter = 100;

cg<decltype(tol)>(I,J,val,x,rhs,N,nz,tol,max_iter);

for(int i=0;i<N;i++)
	{std::cout<< "x[" << i << "]: " << x[i] <<std::endl;}
delete[] x;

delete[] I;
delete[] J;
delete[] val;
delete[] rhs;

std::cout <<"CUDA error: " << cudaGetErrorString(cudaGetLastError()) << std::endl;
return 0;
}
