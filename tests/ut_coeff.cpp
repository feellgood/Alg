#define BOOST_TEST_MODULE coeffTest

#include <boost/test/unit_test.hpp>
#include <vector>
#include <iostream>

#include "alg.h" 

BOOST_AUTO_TEST_SUITE(ut_alg)

BOOST_AUTO_TEST_CASE(affectation_v_coeff)
{
std::cout<< "check implicit operator= on v_coeff" <<std::endl;

alg::v_coeff a(42,3.14);
BOOST_CHECK( a.getVal() == (double) 3.14 );

alg::v_coeff b = a;
BOOST_CHECK( b.getVal() == a.getVal() );

a.setVal(25.0);
BOOST_CHECK( a.getVal() == (double) 25.0 );
BOOST_CHECK( b.getVal() == (double) 3.14 );
}

BOOST_AUTO_TEST_CASE(affectation_by_ref_v_coeff)
{
std::cout<< "check valRef on v_coeff" <<std::endl;

alg::v_coeff a(42,3.14);
a.valRef() = 25.0;
BOOST_CHECK( a.getVal() == (double) 25.0 );

double &val = a.valRef();
val = 4;
BOOST_CHECK( a.getVal() == (double) 4.0 );
alg::v_coeff b = a;
val = 8;
BOOST_CHECK( b.getVal() == (double) 4.0 );
BOOST_CHECK( a.getVal() == (double) 8.0 );
}

BOOST_AUTO_TEST_CASE(operator_comp_v_coeff)
{
std::cout<< "check operator< on v_coeff" <<std::endl;

alg::v_coeff a(42,3.14);
alg::v_coeff b(1,5.0);

BOOST_CHECK( (a < a) == false );
BOOST_CHECK( b < a );
}

BOOST_AUTO_TEST_CASE(equal_index_v_coeff)
{
std::cout<< "check overloaded operator== on v_coeff" <<std::endl;

alg::v_coeff a(42,3.14);
alg::v_coeff b = a;
BOOST_CHECK( b == a );
a.setVal(25.0);
BOOST_CHECK( b == a );
alg::v_coeff c(666,3.14);
BOOST_CHECK( (c == b) == false );// tricky, operator != is not defined for v_coeff
}

BOOST_AUTO_TEST_CASE(affectation_m_coeff)
{
std::cout<< "check implicit operator= on m_coeff" <<std::endl;
alg::m_coeff a(666,42,3.14);
alg::m_coeff b = a;

BOOST_CHECK( a.getVal() == (double) 3.14 );
BOOST_CHECK( b.getVal() == a.getVal() );
a.setVal(25.0);
BOOST_CHECK( a.getVal() == (double) 25.0 );
BOOST_CHECK( b.getVal() == (double) 3.14 );
}

BOOST_AUTO_TEST_CASE(equal_index_m_coeff)
{
std::cout<< "check overloaded operator== on m_coeff" <<std::endl;

alg::m_coeff a(1,42,3.14);
alg::m_coeff b = a;
BOOST_CHECK( b == a );
a.setVal(25.0);
BOOST_CHECK( b == a );
alg::m_coeff c(1,666,3.14);
BOOST_CHECK( (c == b) == false );// tricky, operator != is not defined for m_coeff
}

BOOST_AUTO_TEST_SUITE_END()
