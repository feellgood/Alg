PROJECT (ALG)

find_package (Boost COMPONENTS system filesystem unit_test_framework REQUIRED)

include_directories (${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src
  ${Boost_INCLUDE_DIRS}
  )

set( BOOST_TEST_LOG_LEVEL "all" )

add_definitions (-DBOOST_TEST_DYN_LINK)

add_executable( test_ut_coeff ut_coeff.cpp )

add_executable( test_ut_sparseVect ut_sparseVect.cpp )

add_executable( test_ut_sparseMat ut_sparseMat.cpp )

add_executable( test_ut_core ut_core.cpp )

add_executable( test_ut_cg_dir ut_cg_dir.cpp )

target_link_libraries(test_ut_coeff -ltbb
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

target_link_libraries(test_ut_sparseVect -ltbb
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

target_link_libraries(test_ut_sparseMat -ltbb
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

target_link_libraries(test_ut_core -ltbb
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

target_link_libraries(test_ut_cg_dir -ltbb
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )

add_test (NAME ut_coeff_algebra COMMAND test_ut_coeff)
add_test (NAME ut_sparseVect_algebra COMMAND test_ut_sparseVect)
add_test (NAME ut_sparseMat_algebra COMMAND test_ut_sparseMat)
add_test (NAME ut_core COMMAND test_ut_core)
add_test (NAME ut_cg_dir_basic_pb COMMAND test_ut_cg_dir)
