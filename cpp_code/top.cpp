#include "check.h"
#include "OC.h"
#include "FE.h"
#include "top.h"

using namespace Eigen;
using namespace std;


MatrixXd top(unsigned int nelx, unsigned int nely, double volfrac, double penal, double rmin) {
	int loop = 0; //Used to count the number of iterations in the output. 
	double c = 0;
	double change = 1.0; //Set to the maximum change between x and xold (convergence criteria). 
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
	double length = 20; //! Length
	double breadth = 10; //! Breadth
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
    double ch_change = 0.0;
    double move_ = 0.2;
	while (change > 0.01) {
//    while(loop < 15){
		
		loop++;
		
		xold = x;

		fe_object.assemble(x,penal);
		U = fe_object.solve();
        fe_object.fem_to_vtk();
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
//                        std::cout<<Ue<<std::endl;
						// Issue in line below (1x4) * (8x8) x (4x1) - needs to be same dimension
						mat_res = Ue.transpose() * fe_object.Klocal * Ue;
                		// FE implementation is all in mat_res
                		c += pow(x(ely, elx ), penal)* mat_res; //*(transpose of Ue) * KE * Ue
						dc(ely, elx) = -penal * pow(x(ely, elx), (penal - 1.))*mat_res; //*(transpose of Ue) * KE * Ue;
			}
		}
        fe_object.saveData("ue.csv", Ue);
//        std::cout<<dc<<std::endl;
        fe_object.saveData("dc.csv", dc);
		// Function check below causes Eigen errors (has to be something small in the function)
		dc = check(nelx, nely, rmin, x, dc);
        fe_object.saveData("dcn.csv", dc);
		x = OC(nelx, nely,volfrac, x, dc);
        fe_object.saveData("xx.csv", x);

		//cout << x << endl;

		xchange = (x - xold).cwiseAbs();

		//printf("\n");
//		cout << xchange << endl;

		change = xchange.maxCoeff();
//        if(ch_change == change){
//            move_ = move_/2;
//        }
//        ch_change = change;
        

		printf("\nIteration # %d, change = %f\n",loop,change);
        cout << x << endl;
	}

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
