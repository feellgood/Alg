#ifndef ALG_CG_H
#define ALG_CG_H

#include "alg_core.h"
#include "alg_sparseMat.h"
#include <iostream>


/** if (T == true) Y = Op(B,A*X) else Y = Op(A*X,B) with sparseMat A, Operator Op will act on the result of A*X and B component to component */
template<bool T>
void LinComb(alg::sparseMat const& A,std::vector<double> const& X,std::vector<double> const&B,std::vector<double> &Y,const std::function<double(double,double)> Op)
{
const size_t _size = X.size();
Y.resize(_size);
if (A.getDim() == _size) A.mult(X,Y);

if (T) std::transform(std::execution::par,Y.begin(),Y.end(),B.begin(),Y.begin(),[Op](double _x,double _y){return Op(_x,_y);} );
else std::transform(std::execution::par,Y.begin(),Y.end(),B.begin(),Y.begin(),[Op](double _x,double _y){return Op(_y,_x);} );
}



namespace alg
{

/** conjugate gradient with diagonal preconditioner, returns residu */

double cg(alg::sparseMat& A, std::vector<double> & x,const std::vector<double> & b, alg::iteration &iter) 
{
double rho, rho_1(0.0);
const size_t DIM = x.size();
std::vector<double> p(DIM),q(DIM),r(DIM),z(DIM), diag_precond(DIM);    

// le preconditionneur diagonal est une matrice diagonale contenant les inverses des coeffs de diag(A), ici on va stocker les coefficients dans un std::vector
A.buildDiagPrecond(diag_precond);

iter.set_rhsnorm(alg::norm(b));

LinComb<false>(A,x,b,r,std::minus<double>()); // r = b - A x 

alg::p_direct(diag_precond,r,z);//mult(P, r, z);
rho = alg::dot(z,r);//rho = vect_sp(z, r);
p.assign(z.begin(),z.end());//copy(z, p);

while (!iter.finished_vect(r)) {
      if (!iter.first()) { 
 	        alg::p_direct(diag_precond,r,z);//mult(P, r, z);
	        rho = alg::dot(z,r);//rho = vect_sp(z, r);
	        alg::scaled(rho/rho_1,p); // p *= (rho/rho1)
		alg::add(z,p);// p += z	; donc  p = z + (rho/rho_1)*p        
		}
      alg::mult(A, p, q);
          
	double a=rho/alg::dot(q,p); //a = rho / vect_sp(q, p);	
	alg::scaled_add(p, +a, x); //add(scaled(p, +a), x);
	alg::scaled_add(q, -a, r);//add(scaled(q, -a), r);
      rho_1 = rho;
      ++iter;
          }   
return alg::norm(r)/alg::norm(b);
}

}//end namespace alg

#endif
