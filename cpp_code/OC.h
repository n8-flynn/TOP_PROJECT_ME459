#ifndef OC_h
#define OC_h

#include <iostream> 
#include <algorithm>
#include <cmath>
#include <stddef.h>
#include "Eigen"

using namespace Eigen;

MatrixXd OC(size_t nelx, size_t nely, double volfrac, MatrixXd& x, MatrixXd& dc);

#endif
