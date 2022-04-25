#ifndef ALG_UTILS_H
#define ALG_UTILS_H

namespace alg
{

/** return square */
template <typename T>
inline T sq(T x) {return x*x;}

/** multiplication of a sparse matrix with a dense vector ; CSR sparse matrix  */
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

/**
build diagonal preconditionner of CSR matrix ; all diagonal coefficients must be non zero
diagDim might be different from N
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
void leftMult_diagPrecond_CSRmat(int *I,int *J, T *val, T *x,const int N, T * diag)
{
for(int i=0;i<N;i++)
	for(int k=I[i];k<I[i+1];k++)
		{ val[k] *= diag[i]; }
}

} // end namespace alg

#endif

