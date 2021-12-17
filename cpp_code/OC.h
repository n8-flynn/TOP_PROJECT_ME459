#ifndef OC_h
#define OC_h

#include <iostream> 
#include <algorithm>
#include <cmath>
#include <stddef.h>
#include "Eigen"

using namespace Eigen;
using namespace std;

MatrixXd OC(unsigned int nelx, unsigned int nely, double volfrac, MatrixXd x, MatrixXd dc, double move_val);

MatrixXd mmax(MatrixXd m1, MatrixXd m2, MatrixXd m3);

MatrixXd mmin(MatrixXd m1, MatrixXd m2, MatrixXd m3);

void writeToCsv2(string fileName, MatrixXd  matrix);

#endif
