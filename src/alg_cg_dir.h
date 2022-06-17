#include "alg.h"
#include "alg_sparseMat.h"


/**
 if (T == true) Y = Op(B,A*X) else Y = Op(A*X,B) with sparseMat A, Operator Op will act on the result of A*X and B component to component with respect to mask b */
template<bool T>
void maskedLinComb(std::vector<size_t> const& v_idx, std::vector<bool> const &b, alg::sparseMat const& A,std::vector<double> const& X,std::vector<double> const&B,std::vector<double> &Y,const std::function<double(double,double)> Op)
{
const size_t _size = X.size();
Y.resize(_size);
if (A.getDim() == _size) A.maskedMult(b,X,Y);

if (T) std::transform(std::execution::par,Y.begin(),Y.end(),B.begin(),Y.begin(),[Op](double _x,double _y){return Op(_x,_y);} );
else std::transform(std::execution::par,Y.begin(),Y.end(),B.begin(),Y.begin(),[Op](double _x,double _y){return Op(_y,_x);} );

alg::zeroFill(v_idx,Y);
}

namespace alg
{
/** conjugate gradient with diagonal preconditioner with Dirichlet conditions, returns residu */

double cg_dir(alg::sparseMat& A, std::vector<double> & x,const std::vector<double> & rhs, const std::vector<double> & xd, const std::vector<size_t>& ld, alg::iteration &iter) 
{
double rho, rho_1(0.0);
const size_t DIM = x.size();
if (rhs.size()!=DIM){std::cout << "rhs size mismatch" << std::endl; exit(1);}

std::vector<double> p(DIM),q(DIM),r(DIM),z(DIM),diag_precond(DIM), b(DIM);    
std::vector<bool> mask;
alg::buildMask(DIM,ld,mask);

A.buildDiagPrecond(mask,diag_precond);

//alg::LinComb<false>(A,xd,rhs,b,std::minus<double>());// b = rhs - A xd
//zeroFill(ld, b);
maskedLinComb<false>(ld,mask,A,xd,rhs,b,std::minus<double>());// b = rhs - A xd with respect to mask

iter.set_rhsnorm(alg::norm(b));
	
//alg::LinComb<false>(A,x,b,r,std::minus<double>()); // r = b - A x
//zeroFill(ld,r);
maskedLinComb<false>(ld,mask,A,x,b,r,std::minus<double>()); // r = b - A x with respect to mask

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
      
      alg::maskedMult(mask,A,p,q);
      
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
