#ifndef OC_h
#define OC_h

#include <stddef.h>
#include <Eigen/Dense>

using namespace Eigen;

MatrixXd OC(size_t nelx, size_t nely, double volfrac, MatrixXd* x, MatrixXd* dc);

#endif