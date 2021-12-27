#include "alg.h"
#include "alg_sparseMat.h"

namespace alg
{



/** conjugate gradient with diagonal preconditioner with Dirichlet conditions, returns residu */

double cg_dir(alg::sparseMat& A, std::vector<double> & x,const std::vector<double> & rhs, const std::vector<double> & xd, const std::vector<size_t>& ld, alg::iteration &iter) 
{
double rho, rho_1(0.0);
const size_t DIM = x.size();
if (rhs.size()!=DIM){std::cout << "rhs size mismatch" << std::endl; exit(1);}

std::vector<double> p(DIM),q(DIM),r(DIM),z(DIM),diag_precond(DIM), b(DIM);    
b.assign(rhs.begin(),rhs.end());// b = rhs;

// le preconditionneur diagonal est une matrice diagonale contenant les inverses des coeffs de diag(A), ici on va stocker les coefficients dans un std::vector
A.buildDiagPrecond(diag_precond);

alg::mult(A, xd, z); 
alg::sub(z, b);      // b = b - A xd

zeroFill(ld, b);
zeroFill(ld,diag_precond);

iter.set_rhsnorm(alg::norm(b));
	
alg::LinComb<false>(A,x,b,r,std::minus<double>()); // r = b - A x

zeroFill(ld,r);

alg::p_direct(diag_precond,r,z);// z = diag_precond*r
rho = alg::dot(z,r);// rho = z.r
p.assign(z.begin(),z.end()); // p = z

while (!iter.finished_vect(r)) 
	{
      if (!iter.first())
      		{ 
 	        alg::p_direct(diag_precond,r,z);// z = diag_precond*r
	        rho = alg::dot(z,r);
	        alg::scaled(rho/rho_1,p); // p *= (rho/rho1)
		alg::add(z,p);// p += z	; donc  p = (rho/rho_1)*p + diag_precond*r        
		}
      alg::mult(A, p, q);
      zeroFill(ld,q);  
      double a=rho/alg::dot(q,p);
	alg::scaled_add(p, +a, x); //add(scaled(p, +a), x);
	alg::scaled_add(q, -a, r);//add(scaled(q, -a), r);
      rho_1 = rho;
      ++iter;
          }   
alg::add(xd, x);//x += xd
return alg::norm(r)/alg::norm(b);
}

}//end namespace alg
