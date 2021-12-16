#!/usr/bin/env zsh
path_to_eigen="$PWD/Eigen"
g++ -I $path_to_eigen main.cpp OC.cpp check.cpp FE.cpp top.cpp -Wall -O3 -o trial
