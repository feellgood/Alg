#ifndef ALG_LU_H
#define ALG_LU_H

#include "alg_core.h"

namespace alg
{

/** LU solver */
void lu_solve(alg::sparseMat& LU, const std::vector<double> & b, std::vector<double> & x);

/** ILU algo, in place */
void ilu(alg::sparseMat& A);

}//end namespace alg

#endif
