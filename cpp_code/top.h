#ifndef top_h
#define top_h

#include <iostream> 
#include <fstream>
#include <algorithm>
#include <cmath>
#include <stddef.h>
#include "Eigen"

using namespace Eigen;

MatrixXd top(unsigned int nelx, unsigned int nely, double volfrac, double penal, double rmin);

MatrixXd mabs(MatrixXd m1);

#endif


