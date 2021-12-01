#include <stdio.h>
#include <stddef.h>
#include "top.h"


int main(int argc, char* argv[]) 
{
	size_t nelx = 20;
	size_t nely = 10;
	double volfrac = 0.5;
	double penal = 2.0;
	double rmin = 1.125; 

	top(nelx, nely, volfrac, penal, rmin);
	
	return 0;
}
