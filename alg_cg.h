#include "alg.h"
#include <iostream>

namespace alg
{

double cg(alg::r_sparseMat& A, std::vector<double> & x, std::vector<double> & b, alg::iteration &iter) 
{
double rho, rho_1(0.0);
const size_t DIM = x.size();
std::vector<double> p(DIM),q(DIM),r(DIM),z(DIM), diag_precond(DIM);    

// le preconditionneur diagonal est une matrice diagonale contenant les inverses des coeffs de diag(A), ici on va stocker les coefficients dans un std::vector
for(unsigned int i=0;i<diag_precond.size();i++)
	{ diag_precond[i] = 1.0/A(i,i); }

iter.set_rhsnorm(alg::norm(b));
	
r.assign(b.begin(),b.end());// r = b;
std::vector<double> v_temp(DIM); 
alg::mult(A,x,v_temp);// v_temp = A x;
alg::sub(v_temp,r);// r -= v_temp; donc r = b - A x;

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
