#ifndef ALG_CG_DIR_H
#define ALG_CG_DIR_H

#include "alg_core.h"
#include "alg_iter.h"

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

namespace alg
{
/** conjugate gradient with diagonal preconditioner with Dirichlet conditions, returns residu */

double cg_dir(alg::sparseMat& A, std::vector<double> & x,const std::vector<double> & rhs, const std::vector<double> & xd, const std::vector<size_t>& ld, alg::iteration &iter);

}//end namespace alg

#endif //ALG_CG_DIR_H
