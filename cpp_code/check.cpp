#include <iostream>
#include <algorithm>
#include <cmath>
#include "check.h"

using namespace Eigen; 
using namespace std; 

MatrixXd check(int nelx, int nely, double rmin, MatrixXd x, MatrixXd dc) {
	
	int rmin_f = floor(rmin);
	int i;
	int j;
	int k;
	int l;
	double fac = 0;
	int val;
	double *sum = new double;

	MatrixXd dcn(nely, nelx); //Dcn is a matrix that is nelx by nely.

	dcn.setConstant(0.0); //Initializes the Dcn matrix to all zeros. 

	for (i = 0; i < nelx; i++) {
		for (j = 0; j < nely; j++) {
			*sum = 0.0; 
			for (k = max(i + 1 - rmin_f, 1);k <= min(i + 1 + rmin_f, nelx); k++) {
				for (l = max(j + 1- rmin_f, 1); l <= min(j + 1 + rmin_f, nely); l++) {
					val = pow(i + 1 - k, 2) + pow(j + 1 - l, 2);
					fac = rmin - sqrt(pow(i + 1 - k, 2) + pow(j + 1 - l, 2));
					*sum = *sum + max(0.0, fac);
					dcn(j, i) += max(0.0, fac) * x(l-1, k-1) * dc(l-1, k-1);
				}
			}
			dcn(j, i) = dcn(j, i) / (x(j, i) * (*sum));
		}
	} 
	return dcn;
}
