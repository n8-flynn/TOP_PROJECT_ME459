#ifndef TOP_FUNCTIONS
#define TOP_FUNCTIONS

#include <stddef.h>
#include <Eigen/Dense>

using namespace Eigen;

double check(const size_t nelx, const size_t nely, MatrixXd x, MatrixXd dc);
double OC(const size_t nelx, const size_t nely, MatrixXd x, MatrixXd dc);

#endif
