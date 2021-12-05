#ifndef TOP_h
#define TOP_h

#include <Eigen>
#include <stddef.h>

using namespace Eigen;
using namespace std;

MatrixXd top(unsigned int nelx, unsigned int nely, double volfrac, double penal, double rmin);

#endif

