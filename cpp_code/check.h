#ifndef CHECK_h
#define CHECK_h

#include <stddef.h>
#include <stdlib.h>
#include "Eigen"

using namespace Eigen;

MatrixXd check(int nelx, int nely, double rmin, MatrixXd x, MatrixXd dc);

/*!
	\fn MatrixXd check(int nelx, int nely, double rmin, MatrixXd x, MatrixXd dc)
	\brief Optimality criteria based optimizer. 
	\param nelx The number of elements in the horizontal direction. 
	\param nely The number of elements in the vertical direction.
	\param rmin The filter size devided by the size of the elment. 
	\param x An array of design variables.
	\param dc (not sure about this one).
*/

#endif
