#include "OC.h"
#include <cmath>
#include <algorithm>

using namespace Eigen;
using namespace std;

MatrixXd OC(size_t nelx, size_t nely, double volfrac,MatrixXd& x,MatrixXd& dc)
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
	
	MatrixXd move (nelx, nely);
	MatrixXd xnew;
	MatrixXd newdc(nelx, nely);
	
	move.setConstant(0.2);

	while (l2 - l1 > optimal)
	{
		lmid = 0.5 * (l2 + l1);
        //Commented this out as it was giving compile error
		//xnew = max(x, xnew);

		MatrixXd m1 = x - move; 
		MatrixXd m2 = x + move;
		newdc = -1 * dc / lmid;
//		Commented this out as it was giving a complie error. According to the matlab code, you need to take the square root of the elements inside the eigen matrix however this command takes the sqaure root of the matrix itself.
		//newdc.sqrt();
		
		MatrixXd m3 = x * newdc;
		
        //Commented this out as it was giving compile error
		//volatile double ree = (max(optimal, max(mmax(m1, nelx, nely), min(1.0, mmin(m2, nelx, nely), m3))));
		//xnew.setConstant(max(optimal, max(mmax(m1, nelx, nely), min(1.0, mmin(m2, nelx, nely), m3))));
		if (xnew.sum() - volfrac * nelx * nely > 0)
		{
			l1 = lmid;
		}
		else
			l2 = lmid;
	} 
	return xnew;
}

//For this function, I need to test this out ahead of time to even make sure it works before running it in the code 
//The function should return the maximum value of the matrix if done correctly. 
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


//For this function, I need to test this out ahead of time to even make sure it works before running it in the code 
//The function should return the min value of the matrix if done correctly. 
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
