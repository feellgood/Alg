#include "fem.h"

void chapeaux(Fem &fem)
{
const int TRI = fem.TRI;
const int SEG = fem.SEG;

/********************* SEGMENTS *******************/
for (int s=0; s<SEG; s++){
   Seg &seg = fem.seg[s];
   size_t NBN = Seg::NBN;
   size_t NPI = Seg::NPI;

// NPI 2
   double u[NPI]   = {-1., 1.};
   double pds[NPI] = { 1., 1.};
  
   double detJ = seg.len/2.;
 
   alg::denseMat nod(2, NBN);
   for (size_t ie=0; ie<NBN; ie++){
       size_t i= seg.ind[ie];
        nod(0,ie) = fem.node[i].x;
        nod(1,ie) = fem.node[i].y;
       }

   for (size_t npi=0; npi<NPI; npi++){
       alg::denseMat a(NBN, 1), X(2, 1);
       a[0] = (1.-u[npi])/2.;
       a[1] = (1.+u[npi])/2.;

/* positions des points de gauss */
       mult(nod, a, X);
       seg.x[npi]=X[0];
       seg.y[npi]=X[1];

       for (size_t ie=0; ie<NBN; ie++)
           seg.a[ie][npi] = a[ie];
       seg.weight[npi]   = detJ * pds[npi];
       }
   }

/****************** TRIANGLES *****************/
for (size_t t=0; t<TRI; t++){
    Tri &tri = fem.tri[t];
    const size_t NBN = Tri::NBN;
    const size_t NPI = Tri::NPI;

// NPI 3
   double u[NPI]   = {  1/6.,  2/3.,  1/6.};
   double v[NPI]   = {  1/6.,  1/6.,  2/3.};
   double pds[NPI] = {  1/6.,  1/6.,  1/6.};

   double detJ = tri.surf*2.;

   alg::denseMat nod(2, NBN);
    for (size_t ie=0; ie<NBN; ie++){
        size_t i= tri.ind[ie];
        nod(0,ie) = fem.node[i].x;
        nod(1,ie) = fem.node[i].y;
	}

    for (size_t npi=0; npi<NPI; npi++){
        alg::denseMat a(NBN, 1), X(2, 1);
        alg::denseMat da(NBN, 2), dadu(NBN, 2), J(2, 2), invJ(2, 2);
        a[0] = 1.-u[npi]-v[npi];
        a[1] = u[npi];
        a[2] = v[npi];

/* positions des points de gauss */
        mult(nod, a, X);
        tri.x[npi]=X[0];
        tri.y[npi]=X[1];

        dadu(0, 0)= -1;   dadu(0, 1)= -1;  
        dadu(1, 0)= +1;   dadu(1, 1)=  0;   
        dadu(2, 0)=  0;   dadu(2, 1)= +1; 
	
        alg::mult(nod, dadu, J);

        double detJ=J(0, 0)*J(1, 1) - J(0, 1)*J(1, 0);
        if (fabs(detJ) < EPSILON){
            cerr << "jacobienne singuliere dans le triangle " << t << endl;
            exit(1);}

        //lu_inverse(J);
        invJ(0, 0)=+J(1, 1)/detJ;         invJ(0, 1)=-J(0, 1)/detJ;
        invJ(1, 0)=-J(1, 0)/detJ;         invJ(1, 1)=+J(0, 0)/detJ;

//        std::cout << invJ(0, 0) << " "<< invJ(0,1) << " " << invJ(1,0) << " " << invJ(1, 1) << std::endl;
//        exit(1);

        alg::mult(dadu, invJ, da);

        for (int ie=0; ie<NBN; ie++){
            tri.a[ie][npi]   = a[ie];
            tri.dadx[ie][npi]= da(ie,0);
            tri.dady[ie][npi]= da(ie,1);
            }
        tri.weight[npi]    = detJ * pds[npi];
    }
}

}
