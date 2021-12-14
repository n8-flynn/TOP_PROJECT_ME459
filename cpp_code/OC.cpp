#include "OC.h"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
using namespace Eigen;
using namespace std;

void writeToCsv2(string fileName, MatrixXd  matrix)
{
    //! https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html

    const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");
    ofstream file(fileName);
    if (file.is_open())
    {
        file << matrix.format(CSVFormat);
        file.close();
    }
}

MatrixXd mmin(MatrixXd m1, MatrixXd m2,MatrixXd m3) {

    int r = m1.rows();
    int c = m1.cols();

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            if (m1(i, j) < m2(i, j)) {
//                m1(i, j) = m1(i, j);
                m3(i,j) = m1(i,j);
            }
            else {
//                m1(i, j) = m2(i, j);
                m3(i, j) = m2(i, j);
            }
        }
    }
    return m3;
}


MatrixXd mmax(MatrixXd m1, MatrixXd m2, MatrixXd m3) {
    int r = m1.rows();
    int c = m1.cols();
    
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            if (m1(i, j) > m2(i, j)) {
//                m1(i, j) = m1(i, j);
                m3(i,j) = m1(i,j);
            }
            else {
//                m1(i, j) = m2(i, j);
                m3(i,j) = m2(i,j);
            }
        }
    }
    return m3;
}



MatrixXd OC(unsigned int nelx, unsigned int nely, double volfrac,MatrixXd x,MatrixXd dc)
{
	double l1 = 0.0;
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
		
		printf(" \n");
	
//		xnew = mmax(op1m, mmax(x-move, mmin(op2m, mmin(x+move, mdot(x,msqrt(-dc/lmid))))));
        MatrixXd r1 = x.array()*((-dc/lmid).cwiseSqrt()).array();
        MatrixXd r2 = x + move;
        MatrixXd r3(nely,nelx);
        r3 = mmin(r2,r1,r3);
        MatrixXd r4(nely,nelx);
        r4 = mmin(r3,op2m,r4);
        MatrixXd r5 = x - move;
        MatrixXd r6(nely,nelx);
        r6 = mmax(r5,r4,r6);
        MatrixXd r7(nely,nelx);
        r7 = mmax(op1m,r6,r7);
        xnew = r7;
        
        
		//cout << xnew << endl;
        writeToCsv2("xnew.csv", xnew);
        
		if (xnew.sum() - volfrac * nelx * nely > 0)
		{
			l1 = lmid;
		}
		else
			l2 = lmid;
		
	}
	//cout << xnew << endl;
	return xnew;
}





MatrixXd msqrt(MatrixXd m1) {
	int r = m1.rows();
	int c = m1.cols();
	//cout << m1 << endl;
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			m1(i, j) = pow(m1(i, j),0.5);
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


