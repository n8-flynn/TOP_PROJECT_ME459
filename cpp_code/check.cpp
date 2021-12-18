#include <iostream>
#include <algorithm>
#include <cmath>
#include "check.h"

using namespace Eigen; 
using namespace std; 

MatrixXd check(int nelx, int nely, double rmin, MatrixXd x, MatrixXd dc) {
	/*!
		\brief Mesh-independency filter 
		*This filter is used to look through the elemnts of x and determine what elements are required for the structure and what elements can be discarded. 
		\param nelx The number of elements in the horizontal direction.
		\param nely The number of elements in the vertical direction.
		\param rmin The filter size devided by the size of the elment.
		\param x An array of design variables.
		\param dc Stiffness matrix from FE function. 
	*/

	int rmin_f = floor(rmin); 

	double val = 0.0;
	double fac = 0.0;

	MatrixXd dcn(nely, nelx); //New sensitivity matrix baised off of x and dc. 
	dcn.setConstant(0.0); //Initializes the Dcn matrix to all zeros. 

	for (int i = 0; i < nelx; i++) { //Sweeps through the horizontal elements. 
		for (int j = 0; j < nely; j++) { //Sweeps through the vertical elements. 
			double sum =  0.0;
			for (int k = max(i + 1 - rmin_f, 1);k <= min(i + 1 + rmin_f, nelx); k++) { //Sweeps horizontally through elments that are 2*round(rmin) away from the initial element. 
				for (int l = max(j + 1 - rmin_f, 1); l <= min(j + 1 + rmin_f, nely); l++) { //Sweeps vertically through elements that are 2*roufn(rmin) away from the intial element. 
					fac = rmin - sqrt(pow(i + 1 - k, 2) + pow(j + 1 - l, 2));
					sum += max(0.0, fac); //Ensures that the values going into the sum are positive. 
					dcn(j, i) += max(0.0, fac) * x(l-1, k-1) * dc(l-1, k-1);
				}
			}
			dcn(j, i) = dcn(j, i) / (x(j, i) * sum);
		}
	} 
	return dcn;
}
