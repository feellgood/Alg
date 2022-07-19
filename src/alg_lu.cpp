#include "alg_lu.h"

void alg::lu_solve(alg::sparseMat& LU, const std::vector<double> & b, std::vector<double> & x)
{
const size_t N=LU.getDim();
x.clear();
if (x.size() != N){ x.resize(N, 0.);}
 
// L y =b
std::vector<double> y(N, 0.0);
for (size_t i=0; i<N; i++) { 
    double scal=LU(i).dot(y);
    y[i]=b[i]-scal; // diag part of L matrix always equals 1
    }

// U x = Y
for (size_t i=N; i-- >0;) { // iterate from NOD-1 downto 0
    double scal=LU(i).dot(x);
    if (LU(i, i) == 0) {std::cout << "LU solver aborting" << std::endl; exit(1);}
    x[i]=(y[i]-scal)/LU(i, i);
    }
}


void alg::ilu(alg::sparseMat& A)
{
A.collect();  //  coefficients are sorted in lexicographic order

size_t N=A.getDim();
for (size_t i=1; i<N; i++){ // row 0 unchanged  
    alg::sparseVect &Ai = A(i);
    //double &Aii= Ai.getValRef(i);
    auto it_ii = std::find_if(Ai.begin(),Ai.end(),[i](alg::v_coeff & coeff){return (coeff._i == i); } ); 
    double &Aii = it_ii->valRef();
    
    double Aii_bak=Aii;
    for (auto it_ik=Ai.begin(); it_ik != Ai.end() && it_ik->_i < i; ++it_ik){
        size_t k=it_ik->_i;
        double  Akk = A(k, k);
        double &Aik = it_ik->valRef(); // reference to A(i,k)
        if (Akk == 0) {std::cout << "Akk Nul" << std::endl; continue;}
        if (Aik == 0) {std::cout << "Aik Nul" << std::endl; continue;}
        Aik /= Akk;   

        auto it = std::find_if(Ai.begin(),Ai.end(),[&k](alg::v_coeff coeff){return (coeff._i > k); } );
        if (it==Ai.end()) {std::cout << "not found" << std::endl; break;}

        for (auto it_ij = it; it_ij != Ai.end(); ++it_ij){
            size_t j=it_ij->_i;
            double  &Aij = it_ij->valRef(); // reference to A(i,j)
            double  Akj = A(k,j);

            Aij -= Aik *Akj;    
            }

        }
    if (Aii==0) Aii=Aii_bak;
    
   }
}
