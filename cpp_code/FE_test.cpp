#include "FE.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>

#define NSEC_PER_SEC 1000000000



int main(){
    struct timespec start, end;
	unsigned short int no_quad_points = 3; // For quad, this means 9
    unsigned short int nelx = 20; // Number of elements along x
    unsigned short int nely = 10; // Number of elements along y
    double length = 1.; // Length
    double breadth = 1.; // Breadth
    double penal = 2.;
//    double youngs_mod = 2 * pow(10,11);
    double youngs_mod = 1.;
    double pois_rat = 0.3;
    double force = 1.;
    double g = 0.;
    Eigen::MatrixXd x; // This is the Vol frac matrix
    // Since nely is the number of rows and nelx is the number of columns
    x.resize(nely,nelx);
    x.setConstant(nely,nelx,0.5);
    FE tr(nelx,nely,length,breadth,youngs_mod,pois_rat);
    tr.mesh(no_quad_points);
    tr.init_data_structs();
    tr.define_boundary_condition(force,g);
    tr.cal_k_local();
    clock_gettime(CLOCK_MONOTONIC, &start);
    tr.assemble(x,penal);
    clock_gettime(CLOCK_MONOTONIC, &end);
    time_t elapsed_sec_1 = end.tv_sec - start.tv_sec;
    long elapsed_nsec_1 = end.tv_nsec - start.tv_nsec;

    double elapsed_total_1 = elapsed_sec_1 + (double)elapsed_nsec_1 / (double)NSEC_PER_SEC;

    printf("Time taken for fe_imp %g \n",elapsed_total_1 * 1000);
    Eigen::VectorXd U;
    clock_gettime(CLOCK_MONOTONIC, &start);
    U = tr.solve();
    clock_gettime(CLOCK_MONOTONIC, &end);
    time_t elapsed_sec_2 = end.tv_sec - start.tv_sec;
    long elapsed_nsec_2 = end.tv_nsec - start.tv_nsec;

    double elapsed_total_2 = elapsed_sec_2 + (double)elapsed_nsec_2 / (double)NSEC_PER_SEC;

    printf("Time taken for solution of linear system %g \n",elapsed_total_2 * 1000);
    tr.fem_to_vtk();

}
