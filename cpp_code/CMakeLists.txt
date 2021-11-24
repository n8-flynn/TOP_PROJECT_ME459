#This will be used to create a cmake version of the project. 
cmake_minimum_required(VERSION 3.2)

#Project title
project(top)

#Adding files
add_executable(top  top.cpp OC.cpp main.cpp check.cpp FE.cpp)

#Compiling command below 

#gcc top.cpp OC.cpp main.cpp check.cpp FE.cpp  - Wall -O3 -o top
