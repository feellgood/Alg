#define BOOST_TEST_MODULE sparseVectTest

#include <boost/test/unit_test.hpp>
#include <vector>
#include <iostream>

#include "alg_core.h"

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

BOOST_AUTO_TEST_CASE(sparseVect_by_array_to_constructor)
{
alg::sparseVect v( {alg::v_coeff(1,3.14),alg::v_coeff(3,42.0),alg::v_coeff(0,-1.2345) } );

BOOST_CHECK(v.isSorted() == false);

BOOST_CHECK(v.exist(0) == true);
BOOST_CHECK(v.getVal(0) == (double)-1.2345);
BOOST_CHECK(v.exist(1) == true);
BOOST_CHECK(v.getVal(1) == (double)3.14);
BOOST_CHECK(v.exist(3) == true);
BOOST_CHECK(v.getVal(3) == (double)42.0);
BOOST_CHECK(v.exist(2) == false);
}


BOOST_AUTO_TEST_CASE(sparseVect_implicit_affectation)
{
std::cout << "********** test implicit operator= " << std::endl;
alg::sparseVect v;

v.push_back(1,3.14);
v.push_back(2,2.0);
v.push_back(0,4.0);

alg::sparseVect a = v;

BOOST_CHECK(a.exist(0) == true);
BOOST_CHECK(a.getVal(0) == (double)4.0);
BOOST_CHECK(a.exist(1) == true);
BOOST_CHECK(a.getVal(1) == (double)3.14);
BOOST_CHECK(a.exist(2) == true);
BOOST_CHECK(a.getVal(2) == (double)2.0);
BOOST_CHECK(a.exist(3) == false);

BOOST_CHECK( a.isSorted() == v.isSorted() );
BOOST_CHECK( a.isCollected() == v.isCollected() );
v.kill(1);
BOOST_CHECK(a.exist(1) == true);
BOOST_CHECK(a.getVal(1) == (double)3.14);
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
	v.push_back(2,0.5);
	v.push_back(3,-0.25);

v.sort();
v.collect();// to set Collected flag to true

std::cout << v << std::endl;

double val = alg::dot(v,x);
std::cout << "val (should be 1.57)=" << val << std::endl;
BOOST_TEST(val == (double)1.57);
}

BOOST_AUTO_TEST_CASE(p_direct,* boost::unit_test::tolerance(1e-15))
{
	std::vector<double> x = {1,-2,3.14,4,5};
	std::vector<double> y = {1.0,0.0,0.5,-0.25,0};
	std::vector<double> result;

alg::p_direct(x,y,result);
BOOST_TEST(result[0] == (double)1.0);
BOOST_TEST(result[1] == (double)0.0);
BOOST_TEST(result[2] == (double)1.57);
BOOST_TEST(result[3] == (double)-1.0);
BOOST_TEST(result[4] == (double)0.0);
}

BOOST_AUTO_TEST_CASE(MapSparseVect)
{
alg::MapSparseVect x;
x.push_back(1, 0.125);
x.push_back(2, 2.25);
x.push_back(5, 48);
BOOST_CHECK(x.getVal(0) == 0);
BOOST_CHECK(x.getVal(1) == 0.125);
BOOST_CHECK(x.getVal(2) == 2.25);
BOOST_CHECK(x.getVal(3) == 0);
BOOST_CHECK(x.getVal(4) == 0);
BOOST_CHECK(x.getVal(5) == 48);
BOOST_CHECK(x.getVal(6) == 0);

std::vector<double> y{-7, 10, -1, 9, -11, 0.5};
BOOST_CHECK(x.dot(y) == 23);
}


BOOST_AUTO_TEST_SUITE_END()
