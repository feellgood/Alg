#include "cg.h"

#include<iostream>
#include<iomanip>

#include <stdlib.h>
#include <time.h>

#include<type_traits>

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
inline T sq(T x) {return x*x;}

template <typename T>
void multCSR_MatVect(int *I,int *J, T *val, T *x,const int N, T *y)
{
for(int i=0;i<N;i++)
	{
	y[i] = 0;
	for(int k=I[i];k<I[i+1];k++)
		{ y[i] += val[ k ]*x[J[k]]; }
	}
}

bool test_multCSR(void)
{
const int N = 5;
const int nb_coeff = 9;
int *I, *J ;
I = new int[N];
J = new int[nb_coeff];
I[0]= 0;
I[1]= 2;
I[2]= 4;
I[3]= 7;
I[4]= 9;
J[0]= 0;
J[1]= 1;
J[2]= 1;
J[3]= 2;
J[4]= 0;
J[5]= 3;
J[6]= 4;
J[7]= 2;
J[8]= 4;

double *val, *x, *y;
val = new double[nb_coeff];
val[0]= 1.0;
val[1]= 4.0;
val[2]= 2.0;
val[3]= 3.0;
val[4]= 5.0;
val[5]= 7.0;
val[6]= 8.0;
val[7]= 9.0;
val[8]= 6.0;

x = new double[N];
y = new double[4];

x[0]=1.0;
x[1]=2.0;
x[2]=3.0;
x[3]=2.0;
x[4]=1.0;

multCSR_MatVect<double>(I,J,val,x,N,y);

bool test = (y[0] == 9.0) && (y[1] == 13) && (y[2] == 27) && (y[3] == 33); 

delete [] I;
delete [] J;
delete [] val;
delete [] x;
delete [] y;

return test;
}

template <typename T>
T check_sol(int *I,int *J, T *val, T *x,int N, T *rhs)
{
T result(0);
T *y;

y = new T[N];

multCSR_MatVect<T>(I,J,val,x,N,y);

for(int k =0;k<N;k++) 
	{ result += sq<T>( y[k] - rhs[k]); }
delete [] y;

return sqrt(result);
}

int main(void)
{
if(test_multCSR())
	{ std::cout << "mult mat vect with CSR indices Ok." <<std::endl; }

const int N =5;
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
int max_iter = 100;
int nb_iter;

res = GPU::cg<decltype(tol)>(I,J,val,x,rhs,N,tol,max_iter,nb_iter);
std::cout << "nb iter = " << nb_iter << "; residu = " << res << std::endl;

decltype(tol) verif = check_sol(I,J,val,x,N,rhs);
std::cout << "verif = " << verif <<std::endl;

delete[] x;

delete[] I;
delete[] J;
delete[] val;
delete[] rhs;

std::cout <<"CUDA error: " << cudaGetErrorString(cudaGetLastError()) << std::endl;
return 0;
}
