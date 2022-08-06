#ifndef ALG_CG_ILU_DIR_H
#define ALG_CG_ILU_DIR_H

#include "alg_iter.h"
#include "alg_lu.h"

namespace alg
{

/** conjugate gradient with ILU preconditionner with Dirichlet conditions, returns residu */

double cg_ilu_dir(alg::r_sparseMat& A, std::vector<double> & x, const std::vector<double> & rhs, const std::vector<double> & xd, const std::vector<size_t>& ld, alg::iteration &iter);

}//end namespace alg

#endif
