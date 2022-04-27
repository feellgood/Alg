#ifndef ALG_UTILS_H
#define ALG_UTILS_H

#include<type_traits>
#include <vector>

#include "alg_coeff.h"

namespace alg
{

/** return square */
template <typename T>
inline T sq(T x) {return x*x;}

/** multiplication of a sparse matrix with a dense vector ; CSR sparse matrix ; return : y , must be zero initialized  */
template <typename T>
void multCSR_MatVect(int *I,int *J, T *val, T *x,const int N, T *y)
{
for(int i=0;i<N;i++)
	for(int k=I[i];k<I[i+1];k++)
		{ y[i] += val[ k ]*x[J[k]]; }
}

/**
build diagonal preconditionner of CSR matrix ; all diagonal coefficients must be non zero
diagDim might be different from the dimension of the matrix
*/
template <typename T>
void build_diagPrecond_CSRmat(int *I,int *J, T *val, T * diag,const int diagDim)
{
for(int i=0;i<diagDim;i++)
	{
	for(int k=I[i];k<I[i+1];k++)
		{ 
		if (J[k] == i)
			{
			if (val[k]==0.0) { std::cout << "Error : zero coeff on the diagonal." << std::endl; exit(1); }
			else
				{ diag[i] = 1.0/val[k]; }
			}
		}
	}
}


/** 
in place left multiplication of a CSR sparse matrix by a diagonal matrix stored in diag
 */
template <typename T>
void leftMult_diagPrecond_CSRmat(int *I, T *val,const int N, T * diag,const int diagDim)
{
// diagDim might be inferior to N
if(diagDim>N)
	{ exit(1); }
else
	{
	for(int i=0;i<diagDim;i++)
		for(int k=I[i];k<I[i+1];k++)
			{ val[k] *= diag[i]; }
	}
}

/**
build the CSR sparse matrix from the vector of m_coeff C ordered, outputs = I,J,val,N
dim_C is the dimension of Im(C)
I,J and val are memory allocated, must be deleted otherwise memory leak.
Zero index based.
*/
template<typename T>
	void buildCSR_sparseMat(const std::vector<alg::m_coeff> C, const int dim_C,int *I,int *J,T *val,int &N)
{
int nb_coeff = C.size();
N = dim_C;
I = new int[dim_C+1];
I[N] = nb_coeff;
J = new int[nb_coeff];
val = new T[nb_coeff];

int i=0;
for(size_t k=0;k<nb_coeff;k++)
	{
	val[k] = (T) C[k].getVal();
	J[k]=C[k]._j;
	if((C[k]._i == i )&&(i < dim_C)) 
		{I[i] = k;  i++;}
	}
}

} // end namespace alg

#endif

