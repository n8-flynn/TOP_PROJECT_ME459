#include <iostream>
double OC(size_t nelx, size_t nely, double volfrac,double *x, double * dc) {
//nelx is the number of elements in the x direction.
//nely is the number of elements in the y direction. 
//volfrac is the desired volume fraction.
//x is an array of densities.
//dc is ??
 
int l1 = 0;
int l2 = 100000;

double optimal = 0.0001; 
double move = 0.2; 

while (l2-l1 > optimal) {
	double lmid = 0.5 * (l2 + l1);
	*xnew = max(0.0

} 

return xnew;
}
