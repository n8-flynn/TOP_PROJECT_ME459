#ifndef OC_h
#define OC_h

#include <iostream> 
#include <algorithm>
#include <cmath>
#include <stddef.h>
#include "Eigen"

using namespace Eigen;

MatrixXd OC(unsigned int nelx, unsigned int nely, double volfrac, MatrixXd& x, MatrixXd& dc);

#endif
