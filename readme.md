# alg : sparse Matrix implementation for gradient conjugate algorithm

Sparse matrix implemented using a representation of each coefficient that includes indices.
Very few operators implemented yet, just the minimum required to implement gradient conjugate algorithm.

Compile with cmake

some unit tests are also implemented, not build by default.

To build them, type:

cmake . -DENABLE_UTESTS=ON

make test

## Dependencies : 
* C++ 11 or above, STL
* cmake 
* BOOST (unit tests, ..)
