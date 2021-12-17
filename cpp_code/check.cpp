#include <iostream>
#include <algorithm>
#include <cmath>
#include "check.h"

using namespace Eigen; 
using namespace std; 

MatrixXd check(int nelx, int nely, double rmin, MatrixXd x, MatrixXd dc) {
	/*!
		\brief Mesh-independency filter 
		\param nelx The number of elements in the horizontal direction.
		\param nely The number of elements in the vertical direction.
		\param rmin The filter size devided by the size of the elment.
		\param x An array of design variables.
		\param dc(not sure about this one).
	*/

	int rmin_f = floor(rmin);

	double val = 0.0;
	double fac = 0.0;

	MatrixXd dcn(nely, nelx); 
	dcn.setConstant(0.0); //Initializes the Dcn matrix to all zeros. 

	for (int i = 0; i < nelx; i++) {
		for (int j = 0; j < nely; j++) {
			double sum =  0.0;
			for (int k = max(i + 1 - rmin_f, 1);k <= min(i + 1 + rmin_f, nelx); k++) {
				for (int l = max(j + 1 - rmin_f, 1); l <= min(j + 1 + rmin_f, nely); l++) {
					fac = rmin - sqrt(pow(i + 1 - k, 2) + pow(j + 1 - l, 2));
					sum += max(0.0, fac);
					dcn(j, i) += max(0.0, fac) * x(l-1, k-1) * dc(l-1, k-1);
				}
			}
			dcn(j, i) = dcn(j, i) / (x(j, i) * sum);
		}
	} 
	return dcn;
}
