#define BOOST_TEST_MODULE sparseVectTest

#include <boost/test/unit_test.hpp>

#include "alg.h" 

BOOST_AUTO_TEST_SUITE(ut_alg)

BOOST_AUTO_TEST_CASE(Stupid)
{
float x = 1.0;
BOOST_CHECK(x != 0.0f);
}

BOOST_AUTO_TEST_CASE(sparseVect)
{
alg::sparseVect v;

v.push_back(1,3.14);

BOOST_CHECK(v.isSorted() == false);
}


BOOST_AUTO_TEST_CASE(sparseVect_kill_zero)
{
alg::sparseVect v;
std::cout << "********** test sparseVect.kill_zero " << std::endl;

if(v.isEmpty()) {std::cout << "v is empty"  << std::endl;}
else {std::cout << "v is not empty"  << std::endl;}
BOOST_CHECK(v.isEmpty());

v.push_back(1,0);
v.kill_zero();
std::cout<<"after kill_zero"<<std::endl;
if(v.isEmpty()) {std::cout << "v is empty"  << std::endl;}
else {std::cout << "v is not empty"  << std::endl;}
BOOST_CHECK(v.isEmpty());
std::cout << "********** end test sparseVect.kill_zero " << std::endl;
}

BOOST_AUTO_TEST_CASE(sparseVect_kill)
{
alg::sparseVect v;

std::cout << "********** test sparseVect.kill " << std::endl;

v.push_back(1,0);
v.push_back(1,42.0);
v.kill(1);
std::cout<<"after kill"<<std::endl;
if(v.isEmpty()) {std::cout << "v is empty"  << std::endl;}
else {std::cout << "v is not empty"  << std::endl;}
BOOST_CHECK(v.isEmpty());
std::cout << "********** end test sparseVect.kill " << std::endl;
}

BOOST_AUTO_TEST_CASE(sparseVect_collect)
{
alg::sparseVect v;

std::cout << "********** test sparseVect.collect " << std::endl;

v.push_back(1,1);
v.push_back(3,100);
v.push_back(0,3.14);
v.push_back(1,41.0);
v.push_back(3,1.101);
v.push_back(5,-3.14);

v.collect();
std::cout<<"after collect"<<std::endl;
if(v.isEmpty()) {std::cout << "v is empty"  << std::endl;}
else {std::cout << "v is not empty"  << std::endl;}
BOOST_CHECK(!v.isEmpty());
double val = v.getVal(1);
std::cout << "val (should be 42)=" <<  val << std::endl;
BOOST_CHECK (val == (double)42.0);

val = v.getVal(3);
std::cout << "val (should be 101.101)=" <<  val << std::endl;
BOOST_CHECK (val == (double)101.101);

val = v.getVal(0);
std::cout << "val (should be 3.14)=" <<  val << std::endl;
BOOST_CHECK (val == (double)3.14);

val = v.getVal(2);
std::cout << "val (should be 0)=" <<  val << std::endl;
BOOST_CHECK (val == (double)0);

std::cout << "********** end test sparseVect.collect " << std::endl;
}

BOOST_AUTO_TEST_CASE(p_scal,* boost::unit_test::tolerance(1e-15))
{
	std::vector<double> x = {1,-2,3.14,4,5};
	alg::sparseVect v;

	v.push_back(0,1.0);
	v.push_back(10,101);
	v.push_back(2,0.5);
	v.push_back(3,-0.25);

v.sort();
v.collect();// to set Collected flag to true

std::cout << v << std::endl;

double val = alg::p_scal(v,x);
std::cout << "val (should be 1.57)=" << val << std::endl;
BOOST_TEST(val == (double)1.57);
}

BOOST_AUTO_TEST_CASE(w_sparseMat1)
{
	alg::w_sparseMat m(5);
	m.push_back(1,1,3.14);
	m.push_back(1,1,1000);
	m.push_back(1,0,1);
	m.push_back(3,1,3.14);
	m.push_back(0,0,42);
	m.push_back(2,2,456);
	m.push_back(3,3,798);
	m.push_back(4,4,159);
	m.push_back(4,0,2.0);
	std::cout << "m=" << m << std::endl;

	alg::r_sparseMat r_m(m);
	r_m.print();
	BOOST_CHECK(m.getDim() == (size_t)5);

	std::vector<double> x = {0.5,0,0,0,1};
	std::vector<double> y;
	alg::mult(r_m,x,y);
	for(unsigned int i=0;i<y.size();i++) {std::cout << i <<':'<<y[i] <<" ; ";}
}

BOOST_AUTO_TEST_SUITE_END()
