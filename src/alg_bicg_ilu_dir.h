#ifndef ALG_BICG_ILU_DIR_H
#define ALG_BICG_ILU_DIR_H

#include "alg_iter.h"
#include "alg_lu.h"

namespace alg
{

/** biconjugate gradient with ILU preconditioner with Dirichlet conditions, returns residu */
double bicg_ilu_dir(alg::sparseMat& A, std::vector<double> & x, const std::vector<double> & rhs, const std::vector<double>& xd, const std::vector<size_t>& ld, alg::iteration &iter);

}//end namespace alg

#endif
