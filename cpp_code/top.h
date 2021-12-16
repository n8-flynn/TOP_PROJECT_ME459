#ifndef top_h
#define top_h

#include <iostream> 
#include <fstream>
#include <algorithm>
#include <cmath>
#include <stddef.h>
#include "Eigen"

using namespace Eigen;

MatrixXd top(unsigned int nelx, unsigned int nely, double volfrac, double penal, double rmin,int wh);

/*! 
\fn MatrixXd top(unsigned int nelx, unsigned int nely, double volfrac, double penal, double rmin,int wh)
\brief Topology optimization function. Calls FE, OC, and Check function.
\param nelx Number of elements in the x direction.
\param nely Number of elements in the y direction.
\param volfrac Volume fraction.
\param penal Assigned penality to function (not sure how to explain this)
\param rmin Filter size.
\param wh Loading and support condition (not 100% sure)
*/

#endif


