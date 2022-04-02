#include "cg.h"

#include<iostream>
#include<iomanip>

int main(void)
{
const int N =10;
int nz = (N-2)*3 + 4;
int *I, *J ;
float *val, *x, *rhs;

I = new int[N+1];
J = new int[nz];
val = new float[nz];

genTridiag<float>(I, J, val, N, nz);
//for(int i=0;i<nz;i++)	{std::cout<< "val[" << i << "]: " << val[i] <<std::endl;}

x = new float[N];
rhs = new float[N];

for (int i = 0; i<N; i++)
	{
	x[i] = i/10.f;
	rhs[i] = 1.0;
	std::cout<< "x[" << i << "]: " << x[i] <<std::endl;
	}

float tol = 1e-3;
int max_iter = 100;

cg<float>(I,J,val,x,rhs,N,nz,tol,max_iter);

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
