#This will be used to create a cmake version of the project. 
cmake_minimum_required(VERSION 3.2)

#Project title
project(top)

#Adding files
add_executable(trial main.cpp check.cpp top.cpp OC.cpp FE.cpp check.h OC.h FE.h top.h)

#Compiling command below 

#gcc main.cpp top.cpp OC.cpp check.cpp FE.cpp - Wall -O3 -o main
