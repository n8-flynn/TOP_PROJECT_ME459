#include <iostream>
#include <algorithm>
#include <Eigen/Dense>
#include <numeric>
#include "OC.h"

using namespace Eigen;

MatrixXd OC(size_t nelx, size_t nely, double volfrac,MatrixXd x,MatrixXd dc) 
{
	//nelx is the number of elements in the x direction.
	//nely is the number of elements in the y direction. 
	//volfrac is the desired volume fraction.
	//x is an array of densities.
	//dc is ??
 
	int l1 = 0;
	int l2 = 100000;

	double optimal = 0.0001; 
	MatrixXd move (nelx, nely);
	move.setConstant(0.2);
	double lmid;
	double sum; 

	MatrixXd xnew; 

	while (l2 - l1 > optimal)
	{
		lmid = 0.5 * (l2 + l1);

		xnew = max(optimal, max(x - move, min(1.0, min(x+move,x*sqrt(-dc/lmid)))));
		if (xnew.sum() - volfrac * nelx * nely > 0)
		{
			l1 = lmid;
		}
		else
			l2 = lmid;
	} 
	return xnew;
}
