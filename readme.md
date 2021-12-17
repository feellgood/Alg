# Alg - A set of linear solvers for sparse matrices

Alg is a library, header only, to handle sparse matrices involved in some finite element problems. It is NOT a general purpose library, only the minimal set of elementary algebraic functions is implemented in order to give access to gradient conjugate and bi-conjugate gradient algorithms.

Available algorithms:
conjugate gradient : standard and directional
biconjugate stabilized gradient : standard and directional
LU(for precond matrix)
diagonal preconditioner is implemented in most routines.

All cg and related algorithms are using diagonal pre-conditioner.

### Dependencies :
C++ 17 and STL

The package libtbb-dev has to be installed because STL depends on it for the implementation of parallel execution policies. You have to link with -ltbb (see CMakeLists.txt).
