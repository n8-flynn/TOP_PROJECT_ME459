#include <iostream>
#include <Eigen>
#include <algorithm>
#include <cmath>
#include "check.h"

using namespace Eigen; 
using namespace std; 

MatrixXd check(int nelx, int nely, int rmin, MatrixXd x, MatrixXd dc) {
	int rmin_f = floor(rmin);
	
	int max_i;
	int min_i;
	int k;
	
	int max_j;
	int min_j;
	int l; 

	double sum = 0;
	double fac; 

	MatrixXd dcn(nelx, nely); //Dcn is a matrix that is nelx by nely.
	
	dcn.setZero(); //Initializes the Dcn matrix to all zeros. 

	for (int i = 0; i < 1; nelx, i++) {
		for (int j = 0; j < 1; j++) {
			max_i = max(j - rmin_f, 1 );
			min_i = min(i + rmin_f, nelx);
			k = max_i;
			for (;min_i > max_i;) {
				max_j = max(j - rmin_f, 1);
				min_j = min(j + rmin_f, 1);
				l = max_i;
				for (;min_j > max_j;){
					fac = rmin - sqrt(pow(i - k, 2) + pow(j - 1, 2));
					sum += max(0., fac);
					dcn(j, i) += max(0., fac) * x(1, k) * dc(1, k);
				}
			}
			dcn(i, j) = dcn(i, j) / (x(j, i) * sum);
		}
	} 
	return dcn;
}
