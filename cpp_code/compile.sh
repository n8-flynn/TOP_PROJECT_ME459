#!/usr/bin/env zsh
#SBATCH --job-name=huz_final_project
#SBATCH --partition=wacc
#SBATCH --cpus-per-task=1
#SBATCH --time=0-00:30:00
#SBATCH -o test.out
#SBATCH --mem=8G

path_to_eigen="$PWD/Eigen"
g++ -I $path_to_eigen main.cpp OC.cpp check.cpp FE.cpp top.cpp -Wall -O3 -o run

