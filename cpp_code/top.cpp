#include "check.h"
#include "OC.h"
#include "FE.h"
#include "top.h"

using namespace Eigen;
using namespace std;

MatrixXd top(unsigned int nelx, unsigned int nely, double volfrac, double penal, double rmin,int wh) {
	printf("Top starting\n");
	int loop = 0.0; //Used to count the number of iterations in the output. 
	double change = 1.0; //Set to the maximum change between x and xold (convergence criteria). 
	double move = 0.2; //Used to move the matrix up or down to force optimization convergence.

	MatrixXd x(nely, nelx); //Oringal volume fraction field 
	MatrixXd dc(nely, nelx);  //
	x.setConstant(volfrac); //Sets all elements in matrix x to the volume fraction variable. 
	MatrixXd xold(nely, nelx);
	MatrixXd xchange(nely, nelx);
	xchange.setZero();
	VectorXd U;

	// Defining material properties required for the FEM code
	// This defines the number of quadrature points we use to integrate our finite dimensional weak form in order to compute K elemental - Over here we use the 3 point guass quad
	// rule as this is sufficient to compute the integration exactly
	
	uint8_t no_quad_points = 3; // Domain dimensions
	double length = nely; // Length
	double breadth = nelx; // Breadth
   	double youngs_mod = 1; // Youngs Modulus of the material
    double pois_rat = 0.3; // Poisons ratio
    double force = -1.; // Force acting on the cantilivered beam
    double g = 0.; // The dirichlet boundary condition - For this problem we fix the left edge of the 2D domain

    // Define the FE class object
    FE fe_object(nelx,nely,length,breadth,youngs_mod,pois_rat);
    // For more details about each FE class method, please view documentation
    
    // Generate the mesh - This basically fills up the nodal coordinate and the element connectiviity matrices - It takes no_of_quad_points as we also fill up the values of the quad rule in this function
    fe_object.mesh(no_quad_points);
    fe_object.init_data_structs(); // Initiliaze all the major datastructures to the right size
    fe_object.define_boundary_condition(force,g,wh); // Define the boundary conditions
    fe_object.cal_k_local();
	printf("Solving");

	while (change > 0.0001 && x.sum() > 0.95 * volfrac * nely * nelx) {
		loop++;
		xold = x;
		fe_object.assemble(x, penal); // Generate K global based on x and penalty
		U = fe_object.solve(); // Determine solution vector U
		unsigned short int ele_no;
		vector<unsigned short int> global_nodes;
		VectorXd Ue;
		Ue.resize(fe_object.dofs_per_ele);

		double mat_res = 0;
		for (unsigned short int ely = 0; ely < nely; ely++) {
			for (unsigned short int elx = 0; elx < nelx; elx++) {
				ele_no = elx * nely + ely;
				global_nodes = fe_object.EC[ele_no]; // EC , the elemental connectivity matrix will give us the global nodes corresponding to the element number
				Ue = U(global_nodes); // Extract the elemental U (Ue) from the global U using the global nodes vector
				mat_res = Ue.transpose() * fe_object.Klocal * Ue; // mat_res you could say is almost like the elemental compliance without the relative densities x multiplied. This is done as directly putting this into the below formula gives a compile error for some reason.
				dc(ely, elx) = -penal * pow(x(ely, elx), (penal - 1.)) * mat_res; // dc is nothing but the change in compliance c with respect to the relative densities x. This is nothing but the sensitivity. 
			}
		}

		dc = check(nelx, nely, rmin, x, dc);

		x = OC(nelx, nely, volfrac, x, dc, move); //Optimization criteria function.

		xchange = (x - xold).cwiseAbs(); //Compares the old volume fraction field with the new field. 
		change = xchange.maxCoeff(); //Finds the larges change in the xchange matrix as a single value. 
		printf(".");
	}
	printf("\n Done\n");
	return x;
}
