#ifndef ALG_H
#define ALG_H

/** \file alg.h 
 * \brief set of class to handle sparse matrix operations for gradient conjugate algorithm
 * a sparse vector class
 * a read and a write sparse matrix class
 * various vector operations, scalar and direct products; ...
 * most of the member functions of the classes are using lambdas and C++11 algorithm and numeric
 * */

#include <iostream>
#include <cmath> // sqrt,fabs
#include <algorithm>
#include <numeric> // inner_product

#include "alg_coeff.h"
#include "alg_sparseVect.h"
#include "alg_sparseMat.h"
#include "alg_denseMat.h"
#include "alg_iter.h"

// cg and its variations
#include "alg_cg.h"
#include "alg_cg_dir.h"
#include "alg_cg_ilu.h"
#include "alg_cg_ilu_dir.h"

// bicg and its variations
#include "alg_bicg.h"
#include "alg_bicg_dir.h"
#include "alg_bicg_ilu.h"
#include "alg_bicg_ilu_dir.h"

#endif //ALG_H

