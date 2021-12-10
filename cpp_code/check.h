#ifndef CHECK_h
#define CHECK_h

#include <stddef.h>
#include "Eigen"

using namespace Eigen;

MatrixXd check(int nelx, int nely, int rmin, MatrixXd x, MatrixXd dc);

#endif
