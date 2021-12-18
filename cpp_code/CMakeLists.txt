#This will be used to create a cmake version of the project. 
cmake_minimum_required(VERSION 3.2)

#Project title
project(top)

#Set Cmake build type to by defualt to release
if(NOT CMAKE_BUILD_TYPE)
 set(CMAKE_BUILD_TYPE Release)
endif()

#Set the cmake flags as O3 always
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

#Adding files
add_executable(trial main.cpp top.cpp OC.cpp FE.cpp check.cpp top.h OC.h FE.h check.h)
target_include_directories(trial PRIVATE ${PROJECT_SOURCE_DIR}/Eigen)
#target_link_libraries (trial Eigen3::Eigen)
