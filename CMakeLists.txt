cmake_minimum_required(VERSION 3.0.0)
project(optimization_methods_x64 VERSION 0.1.0)

set(SOURCES ${PROJECT_SOURCE_DIR}/main.cpp)

include(CTest)
enable_testing()
add_subdirectory(unit_tests)

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(/usr/include/eigen3)
add_executable(optimization_methods_x64 main.cpp)

add_library(optimization_methods STATIC ${SOURCES})

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
