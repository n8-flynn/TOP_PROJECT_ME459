#ifndef check_h
#define check_h

#include <stddef.h>
#include <Eigen/Dense>

using namespace Eigen;

inline MatrixXd check(int nelx, int nely, int rmin, MatrixXd x, MatrixXd dc);

#endif
