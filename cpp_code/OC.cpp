#include "OC.h"
#include <cmath>
#include <algorithm>

using namespace Eigen;
using namespace std;

MatrixXd OC(unsigned int nelx, unsigned int nely, double volfrac,MatrixXd x,MatrixXd dc)
{
	double l1 = 0;
	double l2 = 100000;
	double lmid;
	double op1 = 0.001;
	double op2 = 1.0;
	
	MatrixXd move(nely, nelx);
	move.setConstant(0.2);

	MatrixXd xnew(nely, nelx);

	MatrixXd newdc(nely, nelx);
	
	MatrixXd op1m(nely, nelx);
	op1m.setConstant(op1);

	MatrixXd op2m(nely, nelx);
	op2m.setConstant(op2);

	while (l2 - l1 > 0.0001)
	{
		lmid = 0.5 * (l2 + l1);
		//cout << x << endl;
		//printf(" \n");
	
		xnew = mmax(op1m, mmax(x-move, mmin(op2m, mmin(x+move, mdot(x,msqrt(-dc/lmid))))));
		//cout << xnew << endl;

		if (xnew.sum() - volfrac * nelx * nely > 0)
		{
			l1 = lmid;
		}
		else
			l2 = lmid;
		
	} 
	return xnew;
}

MatrixXd mmax(MatrixXd m1, MatrixXd m2) {
	int r = m1.rows();
	int c = m1.cols();
	
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			if (m1(i, j) > m2(i, j)) {
				m1(i, j) = m1(i, j);
			}
			else {
				m1(i, j) = m2(i, j);
			}
		}
	}
	return m1;
}

MatrixXd mmin(MatrixXd m1, MatrixXd m2) {
	int r = m1.rows();
	int c = m1.cols();

	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			if (m1(i, j) < m2(i, j)) {
				m1(i, j) = m1(i, j);
			}
			else {
				m1(i, j) = m2(i, j);
			}
		}
	}
	
	return m1; 
}

MatrixXd msqrt(MatrixXd m1) {
	int r = m1.rows();
	int c = m1.cols();

	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			m1(i, j) = pow(m1(i, j), 0.5);
		}
	}
	return m1;
}

MatrixXd mdot(MatrixXd m1, MatrixXd m2) {
	int r = m1.rows();
	int c = m1.cols();

	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			m1(i, j) = m1(i, j) * m2(i, j);
		}
	}
	return m1;
}


