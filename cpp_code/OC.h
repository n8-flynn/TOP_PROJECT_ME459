#ifndef OC_h
#define OC_h

#include <stddef.h>
#include <Eigen>

using namespace Eigen;

MatrixXd OC(size_t nelx, size_t nely, double volfrac, MatrixXd& x, MatrixXd& dc);

double mmax(MatrixXd x, int nelx, int nely);

double mmin(MatrixXd x, int nelx, int nely);

#endif
