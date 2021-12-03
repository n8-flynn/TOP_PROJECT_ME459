#include <iostream>
#include <algorithm>
#include <Eigen/Dense>
#include <numeric>
#include "OC.h"
#include <cmath>
#include <algorithm>

using namespace Eigen;
using namespace std;

inline MatrixXd OC(size_t nelx, size_t nely, double volfrac,MatrixXd x,MatrixXd dc) 
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
		xnew = max(x, xnew);

		MatrixXd m1 = x - move; 
		MatrixXd m2 = x + move;
		MatrixXd newdc = -dc / lmid;
		
		newdc = newdc.pow(0.5);
		
		MatrixXd m3 = x * newdc;
		
		xnew.setConstant(max(optimal, max(mmax(m1, nelx, nely), min(1.0, mmin(m2, nelx, nely), m3))));
		if (xnew.sum() - volfrac * nelx * nely > 0)
		{
			l1 = lmid;
		}
		else
			l2 = lmid;
	} 
	return xnew;
}

inline double mmax(MatrixXd x,int nelx, int nely)
{
	double maxVal = x(0, 0);

	for (int i = 0; i < nely; i++)
	{
		for (int j = 0; i < nelx; j++)
		{
			if (x(j, i) > maxVal)
			{
				maxVal = x(j, i);
			}
		}
	}
	return maxVal;
}

inline double mmin(MatrixXd x, int nelx, int nely)
{
	double minVal = x(0, 0);

	for (int i = 0; i < nely; i++)
	{
		for (int j = 0; i < nelx; j++)
		{
			if (x(j, i) < minVal)
			{
				minVal = x(j, i);
			}
		}
	}
	return minVal;
}
