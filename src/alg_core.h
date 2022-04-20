#ifndef ALG_CORE_H
#define ALG_CORE_H

/** \namespace alg
 * grab altogether many algebraic functions between various vectors and/or matrices
*/

#include <numeric> // inner_product
#include <execution>
#include "alg_sparseMat.h"
#include "alg_denseMat.h"

namespace alg
{
/** fill with some zeroes vector v with indices vector v_idx, usefull for 'dir' algorithm */
inline void zeroFill(std::vector<size_t> const& v_idx, std::vector<double> &v)
	{ std::for_each(std::execution::par_unseq,v_idx.begin(),v_idx.end(),[&v](const size_t _i){ v[_i] = 0.0; }); }

/**  Y = alpha*X */
inline void scaled(const std::vector<double> & X, const double alpha, std::vector<double> & Y) 
	{ std::transform(std::execution::par,X.begin(),X.end(),Y.begin(),[alpha](const double _x){ return alpha*_x; }); }

/**  Y *= alpha  */
inline void scaled( const double alpha, std::vector<double> & Y) 
	{ std::for_each(std::execution::par,Y.begin(),Y.end(),[alpha](double &_x){ _x *= alpha; }); }

/** 
returns scalar product X.Y 
dev note :
forcing transform_reduce with std::execution::par or std::execution::par_unseq is slower than inner_product, so here transform reduce is used without fixing execution policy.
*/
inline double dot(const std::vector<double> & X,const std::vector<double> & Y)
	{return std::transform_reduce(X.begin(),X.end(),Y.begin(),0.0,std::plus<>(),std::multiplies<>()); } //C++17
//	{ return std::inner_product(X.begin(),X.end(),Y.begin(),0.0); } // C++11

/** direct product : component to component product of the two vectors x and Y, returns Z = X⊗Y ; also known as Hadamard product of two vectors */
inline void p_direct(const std::vector<double> & X,const std::vector<double> & Y,std::vector<double> & Z)
	{
	if (X.size() == Y.size())
		{
		Z.resize(X.size()); 
		std::transform(std::execution::par_unseq,X.begin(),X.end(),Y.begin(),Z.begin(),std::multiplies<>()); 
		}
	else
		{ std::cout<<"vector size mismatch in p_direct.\n"; exit(1);}
	}

// in C++ 17 trans cannot be a template, only available in C++20 and could be template<class std::function<double(double,double)> Op>
/** transformation of Y with operator Op onto X,Y, convenient function to define add, sub scaled_add, or homemade transformations with a lambda. Execution policy is set to parallel  */
inline void trans(std::vector<double> const& X, std::vector<double> & Y, const std::function<double(double,double)> &Op)
{ std::transform(std::execution::par,Y.begin(),Y.end(),X.begin(),Y.begin(),Op  ); }

/** Y += X       */
inline void add(const std::vector<double> & X, std::vector<double> & Y)
	{ trans(X,Y,std::plus<double>() ); }

/** Y -= X       */
inline void sub(const std::vector<double> & X, std::vector<double> & Y)
	{ trans(X,Y,std::minus<double>() ); }

/** Y += alpha*X       */
inline void scaled_add(const std::vector<double> & X,const double alpha, std::vector<double> & Y)
	{ trans(X,Y, [alpha] (const double _x,double _y) { return _x+(alpha*_y); } );}

/** euclidian norm of vector X */
inline double norm(const std::vector<double> & X)
	{ return sqrt(fabs( alg::dot(X,X) )); }

/** Y = A*X with sparseMat A */
inline void mult(alg::sparseMat const& A,std::vector<double> const& X,std::vector<double> &Y)
{
const size_t _size = X.size();
Y.resize(_size);
if (A.getDim() == _size) A.mult(X,Y);
}

/** Y = A*X with sparseMat A and mask b */
inline void maskedMult(std::vector<bool> const &b, alg::sparseMat const& A,std::vector<double> const& X,std::vector<double> &Y)
{
const size_t _size = X.size();
Y.resize(_size);
if (A.getDim() == _size) A.maskedMult(b,X,Y);
}


/** if (T == true) Y = Op(B,A*X) else Y = Op(A*X,B) with sparseMat A, Operator Op will act on the result of A*X and B component to component */
template<bool T>
void LinComb(alg::sparseMat const& A,std::vector<double> const& X,std::vector<double> const&B,std::vector<double> &Y,const std::function<double(double,double)> Op)
{
const size_t _size = X.size();
Y.resize(_size);
if (A.getDim() == _size) A.mult(X,Y);

if (T) std::transform(std::execution::par,Y.begin(),Y.end(),B.begin(),Y.begin(),[Op](double _x,double _y){return Op(_x,_y);} );
else std::transform(std::execution::par,Y.begin(),Y.end(),B.begin(),Y.begin(),[Op](double _x,double _y){return Op(_y,_x);} );
}



/**
 if (T == true) Y = Op(B,A*X) else Y = Op(A*X,B) with sparseMat A, Operator Op will act on the result of A*X and B component to component with respect to mask b */
template<bool T>
void maskedLinComb(std::vector<size_t> const& v_idx, std::vector<bool> const &b, alg::sparseMat const& A,std::vector<double> const& X,std::vector<double> const&B,std::vector<double> &Y,const std::function<double(double,double)> Op)
{
const size_t _size = X.size();
Y.resize(_size);
if (A.getDim() == _size) A.maskedMult(b,X,Y);

if (T) std::transform(std::execution::par,Y.begin(),Y.end(),B.begin(),Y.begin(),[Op](double _x,double _y){return Op(_x,_y);} );
else std::transform(std::execution::par,Y.begin(),Y.end(),B.begin(),Y.begin(),[Op](double _x,double _y){return Op(_y,_x);} );

alg::zeroFill(v_idx,Y);
}

/** Y = A*X with denseMat A */
inline void mult(alg::denseMat const& A,std::vector<double> const& X,std::vector<double> &Y)
{
const size_t _size = X.size();
const size_t ncolsA = A.ncols();
const size_t nrowsA = A.nrows();
assert(ncolsA==_size);

Y.resize(_size);
for (size_t i=0; i<nrowsA; i++) {
       double val(0);
       for (size_t j=0; j<_size; j++) val += A(i,j)*X[j];
       Y[i]=val;
       }
}

/** Y = trans(X)*A */
inline void transposed_mult(std::vector<double> const& X,alg::denseMat const& A,std::vector<double> &Y)
{
const size_t ncolsA = A.ncols();
Y.resize(ncolsA);
const size_t _size = X.size();

for (size_t j=0; j<ncolsA; j++) { 
       double val(0);
       for (size_t i=0; i<_size; i++) val += X[i]*A(i,j);
       Y[j]=val;
       }
}

/** Y = trans(A)*X with denseMat A */
inline void transposed_mult(alg::denseMat const& A,std::vector<double> const& X,std::vector<double> &Y)
{
const size_t _size = X.size();
const size_t ncolsA = A.ncols();
const size_t nrowsA = A.nrows();
assert(nrowsA==_size);

Y.resize(ncolsA);
for (size_t i=0; i<ncolsA; i++) {
       double val(0);
       for (size_t j=0; j<nrowsA; j++) val += A(j, i)*X[j];
       Y[i]=val;
       }
}

/** C = trans(A)*B */
inline void transposed_mult(alg::denseMat const& A,alg::denseMat const& B,alg::denseMat & C)
{
const size_t nrowsA = A.nrows();
const size_t ncolsA = A.ncols();
const size_t ncolsB = B.ncols();
assert(nrowsA==B.nrows());

C.nrows()=ncolsA;
C.ncols()=ncolsB;
C.resize(ncolsA*ncolsB);
for (size_t i=0; i<ncolsA; i++) { 
    for (size_t j=0; j<ncolsB; j++) 
		{
	        double val(0);
	        for (size_t k=0; k<nrowsA; k++) { val += A(k, i)*B(k, j); }
		C(i, j) = val;      
		}
	}
}

/** C = A*B */
inline void mult(alg::denseMat const& A,alg::denseMat const& B,alg::denseMat & C)
{
const size_t nrowsA = A.nrows();
const size_t ncolsA = A.ncols();
const size_t ncolsB = B.ncols();
assert(ncolsA==B.nrows());

C.nrows()=nrowsA;
C.ncols()=ncolsB;
C.resize(nrowsA*ncolsB);
for (size_t i=0; i<nrowsA; i++) { 
    for (size_t j=0; j<ncolsB; j++) 
	{
        double val(0);
        for (size_t k=0; k<ncolsA; k++) { val += A(i, k)*B(k, j); }
	C(i, j) = val;      
	}
     }
}

/** build a logic mask from index list stored in v_idx */
inline void buildMask(const size_t dim,std::vector<size_t> const& v_idx,std::vector<bool> &mask)
	{
	mask.resize(dim);
	std::fill(mask.begin(),mask.end(),true);
	std::for_each(v_idx.begin(),v_idx.end(),[&mask](const size_t _i) {mask[_i] = false;} );
	}

} // end namespace alg

#endif
