#include "OC.h"
#include <cmath>
#include <algorithm>

using namespace Eigen;
using namespace std;

MatrixXd OC(size_t nelx, size_t nely, double volfrac,MatrixXd &x,MatrixXd &dc)
{
	//nelx is the number of elements in the x direction.
	//nely is the number of elements in the y direction. 
	//volfrac is the desired volume fraction.
	//x is an array of densities.
	//dc is ??
 
	double l1 = 0;
	double l2 = 100000;
	double optimal = 0.0001; 
	double lmid;
	
	MatrixXd move(nelx, nely);
	MatrixXd xnew(nelx, nely);
	MatrixXd newdc(nelx, nely);
	
	move.setConstant(0.2);

	while (l2 - l1 > optimal)
	{
		lmid = 0.5 * (l2 + l1);

		xnew.setConstant(max(x.maxCoeff(&nelx, &nely), xnew.maxCoeff(&nelx, &nely)));

		MatrixXd m1 = x - move; 
		MatrixXd m2 = x + move;
		newdc = -1 * dc / lmid;
		
		newdc.array().sqrt();
		
		MatrixXd m3 = x * newdc;
		
		double mmax = m1.maxCoeff(&nelx, &nely);
		double mmin = m2.minCoeff(&nelx, &nely); 
        
		xnew.setConstant(max(optimal, max(mmax, min(1.0, mmin))));
		
		if (xnew.sum() - volfrac * nelx * nely > 0)
		{
			l1 = lmid;
		}
		else
			l2 = lmid;
	} 
	return xnew;
}



