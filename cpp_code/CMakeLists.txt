#This will be used to create a cmake version of the project. 
cmake_minimum_required(VERSION 3.2)

#Project title
project(top)

#Adding files
add_executable(top main.cpp top.cpp top.h OC.cpp OC.h check.cpp check.h FE.cpp FE.h)

#Compiling command below 

#gcc main.cpp top.cpp OC.cpp check.cpp FE.cpp - Wall -O3 -o main
