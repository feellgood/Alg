#include <algorithm>

namespace alg
{
/** 
Y = alpha*X 
*/
inline void scaled(const std::vector<double> & X, const double alpha, std::vector<double> & Y) 
	{ std::transform(X.begin(),X.end(),Y.begin(),[alpha](const double _x){ return alpha*_x; }); }

/** 
Y *= alpha 
*/
inline void scaled( const double alpha, std::vector<double> & Y) 
	{ std::for_each(Y.begin(),Y.end(),[alpha](double &_x){ _x *= alpha; }); }

/**
produit scalaire de X et Y
*/
inline double p_scal(const std::vector<double> & X,const std::vector<double> & Y)
	{ return std::inner_product(X.begin(),X.end(),Y.begin(),0.0); }

/**
produit direct de deux vecteurs : Z = X⊗Y
*/
inline void p_direct(const std::vector<double> & X,const std::vector<double> & Y,std::vector<double> & Z)
	{ for(int i=0;i<Z.size();i++) Z[i]=X[i]*Y[i]; }

/** Y += X       */
inline void inc(const std::vector<double> & X, std::vector<double> & Y)
	{ std::transform(Y.begin(),Y.end(),X.begin(),Y.begin(),std::plus<double>()  ); }

/** Y -= X       */
inline void dec(const std::vector<double> & X, std::vector<double> & Y)
	{ std::transform(Y.begin(),Y.end(),X.begin(),Y.begin(),std::minus<double>()  ); }


/** Y += alpha*X       */
inline void scaled_inc(const std::vector<double> & X,const double alpha, std::vector<double> & Y)
	{ std::transform(Y.begin(),Y.end(),X.begin(),Y.begin(),[alpha] (const double _x,double _y) { return _x+(alpha*_y); }   ); }


/**
norme de X
*/
inline double norm(const std::vector<double> & X)
	{ return sqrt(abs(p_scal(X,X))); }

}
