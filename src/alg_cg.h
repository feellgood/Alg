#ifndef ALG_CG_H
#define ALG_CG_H

#include <iostream>

#include "alg_core.h"
#include "alg_iter.h"


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



namespace alg
{

/** conjugate gradient with diagonal preconditioner, returns residu */

double cg(alg::sparseMat& A, std::vector<double> & x,const std::vector<double> & b, alg::iteration &iter);

}//end namespace alg

#endif
