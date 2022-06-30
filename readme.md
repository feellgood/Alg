# Alg - A set of linear solvers for sparse matrices
Version 1.4.0

Alg is a library to handle sparse matrices involved in some finite element problems. It is NOT a general purpose library, only the minimal set of elementary algebraic functions is implemented in order to give access to gradient conjugate and bi-conjugate gradient algorithms.

The static library libalg.a is created from CMakeLists.txt in alg/src

Available algorithms:
conjugate gradient : standard and directional
biconjugate stabilized gradient : standard and directional
LU(for precond matrix)
diagonal preconditioner is implemented in most routines.

All cg and related algorithms are using diagonal pre-conditioner.

### installation

sudo make install
to copy the static library libalg.a to usr/local/lib and headers to usr/local/include/alg

### Dependencies :
C++ 17 and STL

libtbb-dev has to be installed because STL depends on TBB for the implementation of parallel execution policies with gcc compiler.
how to install :
sudo apt-get install libtbb-dev

### Options : 
1) Some unit tests are available passing the following option to cmake (OFF by default):
cmake . -DENABLE_UTESTS=ON

2) Some GPU executable (dev only) are available using the dedicated option (OFF by default):
cmake . -DENABLE_GPU=ON 
These small executables are here for future development purpose only, they are useless to the alg library at the moment.
