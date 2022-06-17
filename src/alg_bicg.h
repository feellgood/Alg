#ifndef ALG_BICG_H
#define ALG_BICG_H

#include "alg.h"

namespace alg
{

/** biconjugate gradient and diagonal preconditionner */

double bicg(alg::sparseMat& A, std::vector<double> & x, const std::vector<double> & b, alg::iteration &iter);

}//end namespace alg

#endif //ALG_BICG_H
