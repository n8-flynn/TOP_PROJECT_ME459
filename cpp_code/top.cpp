#include "check.h"
#include "OC.h"
#include "FE.h"
#include "top.h"

using namespace Eigen;
using namespace std;

int main()
{
	unsigned int nelx = 20; 
	unsigned int nely = 10; 
	double volfrac = 0.5;
	double penal = 2.0;
	double rmin = 1.125;
	double change = 1.0;
	int loop = 0;
	double c = 0;

	MatrixXd x(nely, nelx);
	MatrixXd dc(nely, nelx);
	x.setConstant(volfrac);
	MatrixXd xold;
	VectorXd U;

	//! Change is the small change in xold and xnew.
	//! Sets the old volume fraction equal to the previous volume fraction x so thast you can compare the two volume fraction.
	//! xold and x are then used to be compared to each other. 
	//! Defining material properties required for the FEM code
	//! This defines the number of quadrature points we use to integrate our finite dimensional weak form in order to compute K elemental - Over here we use the 3 point guass quad
	//! rule as this is sufficient to compute the eintegration exactly
	
	unsigned int no_quad_points = 3; //! Domain dimensions
	double length = 1; //! Length
    double breadth = 1; //! Breadth 
    double youngs_mod = 1; //! Youngs Modulus of the material
    double pois_rat = 0.3; //! Poisons ratio
    double force = 1.; //! Force acting on the cantilivered beam
    double g = 0.; //! The dirichlet boundary condition - For this problem we fix the left edge of the 2D domain

    // Define the FE class object
    FE fe_object(nelx,nely,length,breadth,youngs_mod,pois_rat);
    // Generate the mesh - This basically fills up the nodal coordinate and the element connectiviity matrices - It takes no_of_quad_points as we also fill up the 
    // values of the quad rule in this function
    fe_object.mesh(no_quad_points); //! Initiliaze all the major datastructures to the right size
    fe_object.init_data_structs(); //! Define the boundary conditons - Currently only supports the cantilivered boundary conditions
    fe_object.define_boundary_condition(force,g);
    fe_object.cal_k_local();

	
	while (change > 0.01) {
		loop++;
		xold = x;

		fe_object.assemble(x,penal);
		U = fe_object.solve();
        unsigned int ele_no;
        std::vector<int> global_nodes;
        VectorXd Ue;
        Ue.resize(fe_object.dofs_per_ele);
        double mat_res = 0;
		for (int ely = 0; ely < nely; ely++) {
			for (int elx = 0; elx < nelx; elx++) {
                ele_no = ely * nelx + elx;
                global_nodes = fe_object.EC[ele_no];
                Ue = U(global_nodes);
                mat_res = Ue.transpose()*fe_object.Klocal * Ue;
                // FE implementation is all in mat_res
                c += pow(x(ely, elx), penal)* mat_res; //*(transpose of Ue) * KE * Ue
				dc(ely, elx) = -penal * pow(x(ely, elx), (penal - 1)); //*(transpose of Ue) * KE * Ue;
			}
		}
		dc = check(nelx, nely, rmin, x, dc);
		x = OC(nelx, nely,volfrac, x, dc);
		MatrixXd xchange = (x - xold);
		change = xchange.cwiseAbs().maxCoeff();
	}
}

