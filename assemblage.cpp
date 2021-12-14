#include "fem.h"

/*
void
assemblage(Tri &tri,
           alg::denseMat    &Ke, vector <double> &Le,
           alg::w_sparseMat &K,  vector <double> &L)
{
const size_t NBN = Tri::NBN;

for (size_t ie=0; ie<NBN; ie++){
    size_t i= tri.ind[ie];             
    for (size_t je=0; je<NBN; je++){
        size_t j= tri.ind[je];
        K.push_back(i, j, Ke(ie, je));  
        }
    L[i]+= Le[ie];
    }
}


void
assemblage(Seg &seg,
           alg::denseMat    &Ke, vector <double> &Le,
           alg::w_sparseMat &K,  vector <double> &L)
{
const size_t NBN = Seg::NBN;

for (size_t ie=0; ie<NBN; ie++){
    size_t i= seg.ind[ie];             
    for (size_t je=0; je<NBN; je++){
        size_t j= seg.ind[je];
        K.push_back(i, j, Ke(ie, je));  
        }
    L[i]+= Le[ie];
    }
}
*/
