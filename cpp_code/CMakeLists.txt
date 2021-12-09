#This will be used to create a cmake version of the project. 
cmake_minimum_required(VERSION 3.2)

#Project title
project(top)

#Adding files
add_executable(trial main.cpp top.cpp OC.cpp FE.cpp check.cpp top.h OC.h FE.h check.h)
 
