#include "check.h"
#include "OC.h"
#include "FE.h"
#include "top.h"

using namespace Eigen;
using namespace std;


MatrixXd top(unsigned int nelx, unsigned int nely, double volfrac, double penal, double rmin) {
	printf("Top starting\n");

	int loop = 0.0; //Used to count the number of iterations in the output. 
	double c = 0.0;
	double change = 1.0; //Set to the maximum change between x and xold (convergence criteria). 
	double move = 0.2;

	MatrixXd x(nely, nelx);
	MatrixXd dc(nely, nelx); 
	x.setConstant(volfrac); //Sets all elements in matrix x to the volume fraction variable. 
	MatrixXd xold(nely, nelx);
	MatrixXd xchange(nely, nelx);
	xchange.setZero();
	VectorXd U;

	//! Defining material properties required for the FEM code
	//! This defines the number of quadrature points we use to integrate our finite dimensional weak form in order to compute K elemental - Over here we use the 3 point guass quad
	//! rule as this is sufficient to compute the eintegration exactly
	
	uint8_t no_quad_points = 3; //! Domain dimensions
	double length = nely; //! Length
	double breadth = nelx; //! Breadth
   	double youngs_mod = 1; //! Youngs Modulus of the material
    double pois_rat = 0.3; //! Poisons ratio
    double force = 20.; //! Force acting on the cantilivered beam
    double g = 0.; //! The dirichlet boundary condition - For this problem we fix the left edge of the 2D domain

    // Define the FE class object
    FE fe_object(nelx,nely,length,breadth,youngs_mod,pois_rat);
    // Generate the mesh - This basically fills up the nodal coordinate and the element connectiviity matrices - It takes no_of_quad_points as we also fill up the 
    // values of the quad rule in this function
    fe_object.mesh(no_quad_points); //! Initiliaze all the major datastructures to the right size
    fe_object.init_data_structs(); //! Define the boundary conditons - Currently only supports the cantilivered boundary conditions
    fe_object.define_boundary_condition(force,g);
    fe_object.cal_k_local();
	printf("Solving");

	while (change > 0.0001 && x.sum() > 0.9 * volfrac * nely * nelx) {
		loop++;
		xold = x;
		fe_object.assemble(x, penal);
		U = fe_object.solve();
		unsigned short int ele_no;
		vector<unsigned short int> global_nodes;
		VectorXd Ue;
		Ue.resize(fe_object.dofs_per_ele);

		double mat_res = 0;
		for (unsigned short int ely = 0; ely < nely; ely++) {
			for (unsigned short int elx = 0; elx < nelx; elx++) {
				ele_no = elx * nely + ely;
				global_nodes = fe_object.EC[ele_no];
				Ue = U(global_nodes);
				mat_res = Ue.transpose() * fe_object.Klocal * Ue; // FE implementation is all in mat_res
				c += pow(x(ely, elx), penal) * mat_res; //*(transpose of Ue) * KE * Ue
				dc(ely, elx) = -penal * pow(x(ely, elx), (penal - 1.)) * mat_res; //*(transpose of Ue) * KE * Ue;
			}
		}

		dc = check(nelx, nely, rmin, x, dc);

		x = OC(nelx, nely, volfrac, x, dc, move); //Optimization criteria function.

		xchange = (x - xold).cwiseAbs(); //Compares the old volume fraction field with the new field. 
		change = xchange.maxCoeff(); //Finds the larges change in the xchange matrix as a single value. 
		printf(".");
	}
	printf("\nDone\n");
	return x;
}

MatrixXd mabs(MatrixXd m1) {
	int r = m1.rows();
	int c = m1.cols();

	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			m1(i, j) = pow(pow(m1(i, j), 2),0.5);
		}
	}
	return m1;
}
