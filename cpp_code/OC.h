#ifndef OC_h
#define OC_h

#include <stddef.h>
#include <Eigen/Dense>

using namespace Eigen;

inline MatrixXd OC(size_t nelx, size_t nely, double volfrac, MatrixXd* x, MatrixXd* dc);

inline double mmax(MatrixXd x, int nelx, int nely);

inline double mmin(MatrixXd x, int nelx, int nely);

#endif