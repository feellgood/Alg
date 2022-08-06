#ifndef ALG_CG_ILU_H
#define ALG_CG_ILU_H

#include "alg_iter.h"
#include "alg_lu.h"

namespace alg
{

/** conjugate gradient with ILU preconditioner, returns residu */

double cg_ilu(alg::r_sparseMat& A, std::vector<double> & x, const std::vector<double> & b, alg::iteration &iter);

}//end namespace alg

#endif
