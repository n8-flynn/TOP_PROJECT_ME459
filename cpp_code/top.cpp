#include <iostream> 
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include "check.h"
#include "OC.h"
#include "FE.h"
#include "top.h"

using namespace Eigen;
using namespace std;

inline MatrixXd top(int nelx, int nely, double volfrac, double penal, double rmin) {
	
	MatrixXd x(nelx, nely);
	x.setConstant(volfrac);
	MatrixXd xold;
	MatrixXd dc;
	VectorXd U;

	double change = 1.0;
	int loop = 0;
	//! Change is the small change in xold and xnew.
	//! Sets the old volume fraction equal to the previous volume fraction x so thast you can compare the two volume fraction.
	//! xold and x are then used to be compared to each other. 
	double n1;
	double n2;
	double Ue;
	double c = 0;
	MatrixXd dc(nely, nelx);
	
	while (change > 0.01) {
		loop++;
		xold = x;

		U = FE(nelx, nely)

		for (int ely = 0; ely < nely; ely++) {
			for (int elx = 0; elx < nelx; elx++) {
				
				/* 
				Here goes the FE function. The manipulation of the stiffness matrix goes here. 
				Ue = ?? Up to you Huzaifa. 
				*/
				c += pow(x(ely, elx), penal); //*(transpose of Ue) * KE * Ue
				dc(ely, elx) = -penal * pow(x(ely, elx), (penal - 1)); //*(transpose of Ue) * KE * Ue;
			}
		}
		dc = check(nelx, nely, rmin, x, dc);
		x = OC(nelx, nely,volfrac, &x, &dc);
		MatrixXd xchange = (x - xold);
		change = xchange.cwiseAbs().maxCoeff();
	}
	return x;
}




