#include "OC.h"
#include <cmath>
#include <algorithm>

using namespace Eigen;
using namespace std;

MatrixXd OC(unsigned int nelx, unsigned int nely, double volfrac,MatrixXd &x,MatrixXd &dc)
{
	double l1 = 0;
	double l2 = 100000;
	double lmid;
	
	MatrixXd move(nely, nelx);
	MatrixXd xnew(nely, nelx);
	MatrixXd newdc(nely, nelx);
	MatrixXd m3(nely, nelx);
	
	move.setConstant(0.2);
	
	while (l2 - l1 > 0.0001)
	{
		lmid = 0.5 * (l2 + l1); //eq 1

		MatrixXd m1 = x - move; 
		//cout << m1 << endl; //eq 2
		MatrixXd m2 = x + move;
		//cout << m2 << endl; //eq 3
		//cout << dc << endl; //eq 3
		newdc = -dc/lmid;
		//cout << newdc << endl;
		//newdc.sqrt(); //Need to write 
		m3 = x.cwiseProduct(newdc);
		//cout << m3 << endl;//eq 4
		
		double mmax = m1.maxCoeff(&nely, &nelx); //eq 5
		double mmin = m2.minCoeff(&nely, &nelx); //eq 6
		double sqrtMin = m3.minCoeff(&nely, &nelx); //eq 7

        xnew.setConstant(max(0.001, max(mmax, min(1.0, min( mmin, sqrtMin)))));
		//cout << xnew << endl; //eq 8
		if (xnew.sum() - volfrac * nelx * nely > 0)
		{
			l1 = lmid;
		}
		else
			l2 = lmid;
		
	} 
	return xnew;
}



