PROJECT (ALG)

find_package (Boost COMPONENTS system filesystem unit_test_framework REQUIRED)

include_directories (${CMAKE_SOURCE_DIR}
  ${Boost_INCLUDE_DIRS}
  )

set( BOOST_TEST_LOG_LEVEL "all" )

add_definitions (-DBOOST_TEST_DYN_LINK)

add_executable( test_ut_sparseVect ut_sparseVect.cpp )

target_link_libraries(test_ut_sparseVect
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

add_test (NAME ut_sparseVect_algebra COMMAND test_ut_sparseVect)
