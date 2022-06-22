#include "alg.h"

namespace alg
{

/** biconjugate gradient with dirichlet condtions and diagonal preconditionner */

double bicg_dir(alg::sparseMat& A, std::vector<double> & x, const std::vector<double> & rhs, const std::vector<double>& xd, const std::vector<size_t>& ld, alg::iteration &iter);

}//end namespace alg
