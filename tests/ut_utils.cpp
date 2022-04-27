#define BOOST_TEST_MODULE utilsTest

#include <boost/test/unit_test.hpp>
#include <iostream>

#include "alg_utils.h"

BOOST_AUTO_TEST_SUITE(ut_utils)

BOOST_AUTO_TEST_CASE(multCSR_MV,* boost::unit_test::tolerance(1e-15))
{
const int N = 5;
const int nb_coeff = 9;
int *I, *J ;
I = new int[N];
J = new int[nb_coeff];
I[0]= 0;
I[1]= 2;
I[2]= 4;
I[3]= 7;
I[4]= 9;
J[0]= 0;
J[1]= 1;
J[2]= 1;
J[3]= 2;
J[4]= 0;
J[5]= 3;
J[6]= 4;
J[7]= 2;
J[8]= 4;

double *val, *x, *y;
val = new double[nb_coeff];
val[0]= 1.0;
val[1]= 4.0;
val[2]= 2.0;
val[3]= 3.0;
val[4]= 5.0;
val[5]= 7.0;
val[6]= 8.0;
val[7]= 9.0;
val[8]= 6.0;

x = new double[N];
y = new double[4]();// initialization to zero

x[0]=1.0;
x[1]=2.0;
x[2]=3.0;
x[3]=2.0;
x[4]=1.0;

alg::multCSR_MatVect<double>(I,J,val,x,N,y);
delete [] I;
delete [] J;
delete [] val;
delete [] x;

bool test = (y[0] == 9.0) && (y[1] == 13) && (y[2] == 27) && (y[3] == 33); 
if (test) std::cout << "mult mat vect with CSR indices Ok." <<std::endl;

BOOST_TEST(y[0] == (double)9.0);
BOOST_TEST(y[1] == (double)13.0);
BOOST_TEST(y[2] == (double)27.0);
BOOST_TEST(y[3] == (double)33.0);

delete [] y;
}

BOOST_AUTO_TEST_CASE(diagPrecond,* boost::unit_test::tolerance(1e-15))
{
const int N = 5;
const int nb_coeff = 11;
int *I, *J ;
I = new int[N];
J = new int[nb_coeff];
I[0]= 0;
I[1]= 2;
I[2]= 4;
I[3]= 8;
I[4]= 11;

J[0]= 0;
J[1]= 1;
J[2]= 1;
J[3]= 2;
J[4]= 0;
J[5]= 2;
J[6]= 3;
J[7]= 4;
J[8]= 2;
J[9]= 3;
J[10]= 4;

double *val, *diag;
val = new double[nb_coeff];
val[0]= 1.0;
val[1]= 4.0;
val[2]= 2.0;
val[3]= 3.0;
val[4]= 5.0;
val[5]= -1.0;
val[6]= 7.0;
val[7]= 8.0;
val[8]= 9.0;
val[9]= -2.0;
val[10]= 6.0;

const int diagDim = N-1;
diag = new double[diagDim];

for (int i=0;i<diagDim;i++) { diag[i]=0; }

alg::build_diagPrecond_CSRmat<double>(I,J,val,diag,diagDim);

for (int i=0;i<diagDim;i++) { std::cout << "diag["<< i << "]= " << diag[i] << std::endl; }

delete [] I;
delete [] J;
delete [] val;

bool test = (diag[0] == 1.0) && (diag[1] == 0.5) && (diag[2] == -1.0) && (diag[3] == -0.5); 
if (test) std::cout << "diag from CSR matrix Ok." <<std::endl;

BOOST_TEST(diag[0] == (double)1.0);
BOOST_TEST(diag[1] == (double)0.5);
BOOST_TEST(diag[2] == (double)-1.0);
BOOST_TEST(diag[3] == (double)-0.5);

delete [] diag;
}

BOOST_AUTO_TEST_CASE(multCSRmat_diagPrecond,* boost::unit_test::tolerance(1e-15))
{
const int N = 5;
const int nb_coeff = 11;
int *I, *J ;
I = new int[N];
J = new int[nb_coeff];
I[0]= 0;
I[1]= 2;
I[2]= 4;
I[3]= 8;
I[4]= 11;

J[0]= 0;
J[1]= 1;
J[2]= 1;
J[3]= 2;
J[4]= 0;
J[5]= 2;
J[6]= 3;
J[7]= 4;
J[8]= 2;
J[9]= 3;
J[10]= 4;

double *val, *diag;
val = new double[nb_coeff];
val[0]= 1.0;
val[1]= 4.0;
val[2]= 2.0;
val[3]= 3.0;
val[4]= 5.0;
val[5]= -1.0;
val[6]= 7.0;
val[7]= 8.0;
val[8]= 9.0;
val[9]= -2.0;
val[10]= 6.0;

const int diagDim = N-1;
diag = new double[diagDim]();

alg::build_diagPrecond_CSRmat<double>(I,J,val,diag,diagDim);

alg::leftMult_diagPrecond_CSRmat<double>(I,val,N,diag,diagDim);

alg::build_diagPrecond_CSRmat<double>(I,J,val,diag,diagDim);// diag precond should be unit matrix after the left multiplication above

bool test = (diag[0] == 1.0) && (diag[1] == 1.0) && (diag[2] == 1.0) && (diag[3] == 1.0); 
if (test) std::cout << "diag is unit matrix." <<std::endl;

BOOST_TEST(diag[0] == (double)1.0);
BOOST_TEST(diag[1] == (double)1.0);
BOOST_TEST(diag[2] == (double)1.0);
BOOST_TEST(diag[3] == (double)1.0);

BOOST_TEST(val[0] == (double)1.0);
BOOST_TEST(val[1] == (double)4.0);
BOOST_TEST(val[2] == (double)1.0);
BOOST_TEST(val[3] == (double)1.5);
BOOST_TEST(val[4] == (double)-5.0);
BOOST_TEST(val[5] == (double)1.0);
BOOST_TEST(val[6] == (double)-7.0);
BOOST_TEST(val[7] == (double)-8.0);
BOOST_TEST(val[8] == (double)-4.5);
BOOST_TEST(val[9] == (double)1.0);
BOOST_TEST(val[10] == (double)-3.0);

delete [] I;
delete [] J;
delete [] val;
delete [] diag;
}

BOOST_AUTO_TEST_CASE(buildCSRmat_fromCOO,* boost::unit_test::tolerance(1e-15))
{
int N(0);
int *I, *J;
double *val;

std::vector<alg::m_coeff> C;// there are 10 coefficients
C.push_back(alg::m_coeff(0,0,1.0));
C.push_back(alg::m_coeff(0,1,4.0));
C.push_back(alg::m_coeff(1,1,2.0));
C.push_back(alg::m_coeff(1,2,3.0));
C.push_back(alg::m_coeff(2,0,5.0));
C.push_back(alg::m_coeff(2,2,-1.0));
C.push_back(alg::m_coeff(2,3,7.0));
C.push_back(alg::m_coeff(2,4,8.0));
C.push_back(alg::m_coeff(3,2,9.0));
C.push_back(alg::m_coeff(3,3,-2.0));
C.push_back(alg::m_coeff(3,4,6.0));

const int dim_C = 4;
size_t nb_coeff = C.size();
I = new int[dim_C+1];
J = new int[nb_coeff];
val = new double[nb_coeff];

alg::buildCSR_sparseMat<double>(C,dim_C,I,J,val,N);

BOOST_TEST(I[0]== (int)0);
BOOST_TEST(I[1]== (int)2);
BOOST_TEST(I[2]== (int)4);
BOOST_TEST(I[3]== (int)8);
BOOST_TEST(I[4]== (int)11);//nb_coeff

BOOST_TEST(J[0]== (int)0);
BOOST_TEST(J[1]== (int)1);
BOOST_TEST(J[2]== (int)1);
BOOST_TEST(J[3]== (int)2);
BOOST_TEST(J[4]== (int)0);
BOOST_TEST(J[5]== (int)2);
BOOST_TEST(J[6]== (int)3);
BOOST_TEST(J[7]== (int)4);
BOOST_TEST(J[8]== (int)2);
BOOST_TEST(J[9]== (int)3);
BOOST_TEST(J[10]== (int)4);

BOOST_TEST(val[0] == (double)1.0);
BOOST_TEST(val[1] == (double)4.0);
BOOST_TEST(val[2] == (double)2.0);
BOOST_TEST(val[3] == (double)3.0);
BOOST_TEST(val[4] == (double)5.0);
BOOST_TEST(val[5] == (double)-1.0);
BOOST_TEST(val[6] == (double)7.0);
BOOST_TEST(val[7] == (double)8.0);
BOOST_TEST(val[8] == (double)9.0);
BOOST_TEST(val[9] == (double)-2.0);
BOOST_TEST(val[10] == (double)6.0);

BOOST_TEST(N == dim_C);

std::cout << "CSR matrix Ok." <<std::endl;

delete [] I;
delete [] J;
delete [] val;
}

BOOST_AUTO_TEST_SUITE_END()
