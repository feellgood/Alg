#ifndef ALG_UTILS_H
#define ALG_UTILS_H

#include<type_traits>
#include <vector>

#include "alg_coeff.h"

namespace alg
{

enum precondType {NONE, DIAG};

template <typename T>
class CSR_mat
	{
	public:
		CSR_mat(const int _N, const int nb_coeff) :N(_N) 
			{
			I = new int[N+1];
			I[N] = nb_coeff;
			J = new int[nb_coeff];
			val = new T[nb_coeff];
			}
		
		~CSR_mat()
			{
			delete[] I;
			delete[] J;
			delete[] val;
			}
		
		int *I;
		int *J;
		T *val;
		const int N;// number of line
	};

/** return square(x) */
template <typename T>
inline T sq(T x) {return x*x;}

/** multiplication of a sparse matrix with a dense vector ; CSR sparse matrix ; return : y , y must be zero initialized  */
template <typename T>
void multCSR_MatVect(CSR_mat<T> const& A,T *x,T *y)
{
for(int i=0;i<A.N;i++)
	for(int k=A.I[i];k<A.I[i+1];k++)
		{ y[i] += A.val[ k ]*x[A.J[k]]; }
}

/**
build diagonal preconditionner of CSR matrix ; all diagonal coefficients must be non zero. If a zero is on the diagonal diag is left unmodified.
diagDim might be different from the dimension of the matrix
*/
template <typename T>
void build_diagPrecond_CSRmat(CSR_mat<T> const& A, T * diag,const int diagDim)
{
for(int i=0;i<diagDim;i++)
	for(int k=A.I[i];k<A.I[i+1];k++)
		{ 
		if ((A.J[k] == i)&&(A.val[k]!=0.0))
			{ diag[i] = 1.0/A.val[k]; }
		}
}


/** 
in place left multiplication of a CSR sparse matrix by a diagonal matrix stored in diag
 */
template <typename T>
void leftMult_diagPrecond_CSRmat(CSR_mat<T> & A, T * diag,const int diagDim)
{
// diagDim might be inferior to N
if(diagDim>A.N)
	{ exit(1); }
else
	{
	for(int i=0;i<diagDim;i++)
		for(int k=A.I[i];k<A.I[i+1];k++)
			{ A.val[k] *= diag[i]; }
	}
}

/**
build the CSR sparse matrix from the vector of m_coeff C ordered, outputs = I,J,val,N
dim_C is the dimension of the ending space of the application C = number of lines
I,J and val are memory allocated, must be deleted otherwise memory leak.
Zero index based.
*/
template<typename T>
void buildCSR_sparseMat(const std::vector<alg::m_coeff> C, const size_t dim_C,int *I,int *J,T *val,int &N)
{
size_t nb_coeff = C.size();
N = dim_C;
I[N] = nb_coeff;

size_t i(0);
for(size_t k=0;k<nb_coeff;k++)
	{
	val[k] = (T) C[k].getVal();
	J[k]=C[k]._j;
	if((C[k]._i == i )&&(i < dim_C)) 
		{
		I[i] = k;  
		i++;
		}
	}
}

} // end namespace alg

#endif

