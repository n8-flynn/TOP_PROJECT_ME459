#include "FE.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>



int main(){
	unsigned int no_quad_points = 3; // For quad, this means 9
    unsigned int nelx = 3; // Number of elements along x
    unsigned int nely = 3; // Number of elements along y
    double length = 1; // Length
    double breadth = 1; // Breadth
    double penal = 2;
    double youngs_mod = 1;
    double pois_rat = 0.3;
    FE tr(nelx,nely,length,breadth,penal,youngs_mod,pois_rat);
    tr.mesh(no_quad_points);
}
