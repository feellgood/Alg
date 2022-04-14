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

} // end namespace alg

#endif

