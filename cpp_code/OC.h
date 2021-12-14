#ifndef OC_h
#define OC_h

#include <iostream> 
#include <algorithm>
#include <cmath>
#include <stddef.h>
#include "Eigen"

using namespace Eigen;

MatrixXd OC(unsigned int nelx, unsigned int nely, double volfrac, MatrixXd x, MatrixXd dc, double move_val);

MatrixXd mmax(MatrixXd m1, MatrixXd m2);

MatrixXd mmin(MatrixXd m1, MatrixXd m2);

MatrixXd msqrt(MatrixXd m1);

MatrixXd mdot(MatrixXd m1, MatrixXd m2);

#endif
