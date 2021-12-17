#ifndef CHECK_h
#define CHECK_h

#include <stddef.h>
#include <stdlib.h>
#include "Eigen"

using namespace Eigen;

MatrixXd check(int nelx, int nely, double rmin, MatrixXd x, MatrixXd dc);

#endif
