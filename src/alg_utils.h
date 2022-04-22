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
*/
template <typename T>
void build_diagPrecond_CSRmat(int *I,int *J, T *val, T *x,const int N, T * diag)
{
for(int i=0;i<N;i++)
	{
	for(int k=I[i];k<I[i+1];k++)
		{ 
		if ((J[k] == i)&&(val[k]!=0.0))
			{ diag[k] = 1.0/val[k]; }
		else
			{ exit(1); }
		}
	}
}


/** in place left multiplication of a CSR sparse matrix by its diagonal preconditionner  ; CSR sparse matrix must not have any zero on the diagonal */
template <typename T>
void leftMult_diagPrecond_CSRmat(int *I,int *J, T *val, T *x,const int N, T * diag)
{
for(int i=0;i<N;i++)
	for(int k=I[i];k<I[i+1];k++)
		{ val[k] *= diag[i]; }
}

} // end namespace alg

#endif

