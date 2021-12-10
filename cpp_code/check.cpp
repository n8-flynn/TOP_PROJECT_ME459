#include <iostream>
#include <algorithm>
#include <cmath>
#include "check.h"

using namespace Eigen; 
using namespace std; 

MatrixXd check(int nelx, int nely, int rmin, MatrixXd x, MatrixXd dc) 
{
	int rmin_f = floor(rmin);
	int max_i;
	int min_i;
	int k;
	int max_j;
	int min_j;
	int l; 
	double fac; 

	MatrixXd dcn(nely, nelx); //Dcn is a matrix that is nelx by nely.
	
	dcn.setZero(); //Initializes the Dcn matrix to all zeros. 

	for (int i = 0; i < nelx; i++) // i increases when it less than nelx
	{
		for (int j = 0; j < nely; j++) 
		{
			double sum = 0.0;

			max_i = max(i - rmin_f, 1 );
			min_i = min(i + rmin_f, nelx);

			k = max_i;
			
			for (int k = max_i;k < min_i; k++) 
			{
				max_j = max(j - rmin_f, 1);
				min_j = min(j + rmin_f, nely);

				for (int l = max_j; l < min_j; l++)
				{
					fac = rmin - sqrt(pow(i - k, 2) + pow(j - l, 2));
					sum += max(0., fac);
					dcn(j, i) += max(0., fac) * x(l , k) * dc(l, k);
				}
			}
			dcn(j, i) = dcn(j, i) / (x(j, i) * sum);
		}
	} 
	return dcn;
}
