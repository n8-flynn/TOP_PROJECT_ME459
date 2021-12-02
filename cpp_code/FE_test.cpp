#include "FE.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>



int main(){
	unsigned int no_quad_points = 3; // For quad, this means 9
    unsigned int nelx = 20; // Number of elements along x
    unsigned int nely = 10; // Number of elements along y
    double length = 1; // Length
    double breadth = 1; // Breadth
    double penal = 2;
//    double youngs_mod = 2 * pow(10,11);
    double youngs_mod = 1;
    double pois_rat = 0.3;
    double force = 1.;
    double g = 0.;
    Eigen::MatrixXd x; // This is the Vol frac matrix
    x.resize(nelx,nely);
    x.setZero(nelx,nely);
    FE tr(nelx,nely,length,breadth,penal,youngs_mod,pois_rat);
    tr.mesh(no_quad_points);
    tr.init_data_structs();
    tr.define_boundary_condition(force,g);
    tr.fe_impl(x);
    Eigen::VectorXd U;
    U = tr.solve();

}
