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
y = new double[4];

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

BOOST_AUTO_TEST_SUITE_END()
