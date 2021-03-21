find_package(Boost REQUIRED COMPONENTS unit_test_framework)

include_directories(${PROJECT_SOURCE_DIR}/include)

set(TEST_SOURCES ${PROJECT_SOURCE_DIR}/unit_tests)
set(TEST_LIBS optimization_methods)

include_directories(/usr/include/eigen3)

add_executable(utilities_test ${TEST_SOURCES}/om_utilities_t)
add_executable(golden_section_test ${TEST_SOURCES}/om_golden_section_t)
add_executable(fibonacci_test ${TEST_SOURCES}/om_fibonacci_t)
add_executable(brent_test ${TEST_SOURCES}/om_brent_t)
add_executable(powell_test ${TEST_SOURCES}/om_powell_t)

target_link_libraries(utilities_test ${TEST_LIBS} ${Boost_LIBRARIES})
target_link_libraries(golden_section_test ${TEST_LIBS} ${Boost_LIBRARIES})
target_link_libraries(fibonacci_test ${TEST_LIBS} ${Boost_LIBRARIES})
target_link_libraries(brent_test ${TEST_LIBS} ${Boost_LIBRARIES})
target_link_libraries(powell_test ${TEST_LIBS} ${Boost_LIBRARIES})

add_test(om_utilities_t utilities_test COMMAND utilities_test)
add_test(om_golden_section_t golden_section_test COMMAND golden_section_test)
add_test(om_fibonacci_t fibonacci_test COMMAND fibonacci_test)
add_test(om_brent_t brent_test COMMAND brent_test)
add_test(om_powell_t powell_test COMMAND powell_test)





