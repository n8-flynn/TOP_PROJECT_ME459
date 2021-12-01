// Created by Huzaifa Mustafa Unjhawala 

#include "FE.h"





// Constructor

FE::FE(unsigned int nelx, unsigned int nely, unsigned int length, unsigned int breadth,double penal, double youngs_mod, double pois_rat){
	L = length;
	B = breadth;
	nelx_ = nelx;
	nely_ = nely;
	E = youngs_mod;
	nu = pois_rat;
	penal_ = penal;
}

// Basis function - Internal function needed for fe implementation
double FE::basis_function(unsigned int node , double xi, double eta){
//Kind of hard coded for a bilinear shape function for wuad and linear for triangular - Based on node number, a formual will be choosen using switch-case, then based on the quad point , xi and eta will be substituted to return the value
    double output;
    switch(node) {
        case 0:
            output = 0.25 * (1-xi) * (1-eta);
            break;
        case 1:
            output = 0.25 * (1-xi) * (1+eta);
            break;
        case 2:
            output = 0.25 * (1+xi) * (1+eta);
            break;
        case 3:
            output = 0.25 * (1+xi) * (1-eta);
            break;
        default:
            std::cout<<"There is no "<<node<<" node, invalid input"<<std::endl;
            break;
    }
	return output;
}
// Basis function defined using the general formula obtained from
std::vector<double> FE::basis_gradient(unsigned int node,double xi, double eta){
// Hard coding the basis gradient as could not derive/find the general formula in the case of 2D - maybe its just a multiplication.
    
    std::vector<double> bg(dim,0.0);
    switch(node) {
        case 0:
            bg[0] = -0.25 * (1 - eta);
            bg[1] = -0.25 * (1 - xi);
            break;
        case 1:
            bg[0] = 0.25 * (1 - eta);
            bg[1] = -0.25 * (1 - xi);
            break;
        case 2:
            bg[0] = 0.25 * (1 - eta);
            bg[1] = 0.25 * (1 - xi);
            break;
        case 3:
            bg[0] = -0.25 * (1 - eta);
            bg[1] = 0.25 * (1 - xi);
            break;
        default:
            std::cout<<"There is no "<<node<<" node, invalid input"<<std::endl;
    }

    return bg;
}

void FE::mesh(unsigned int no_quad_points){
	std::cout<<"Generating Mesh .."<<std::endl;

    // The number of nodes is just an extension of 1D
    
    no_of_nodes = (nelx_ + 1) * (nely_ + 1);
	std::cout<<"Total no. of nodes is "<<no_of_nodes<<std::endl;
    // Each node has 2 degrees of freedom as the number of dimensions is 2
    dim = 2;
    total_dofs = no_of_nodes * dim;

	// Nodal Coordinate array will give the coordinates of each dof not just each node
	NC.resize(total_dofs);
    // Since NC is a vector of a vector , we need to initialize each row of EC
    for(int i = 0; i < total_dofs; i++){
        NC[i] = std::vector<double>(dim,0.0);
    }



	// Make NC - NC remains the same for both the quad and the triangular elements

	double incr_x = L/(nelx_); // since the nodes are equally spaced, the x coordinate will differ by this increment
    double incr_y = L/(nely_); // similarly, the y coordinate will differ by this incremenet
    double x = 0.0; // first node is 0,0
    double y = 0.0;
    // Construct NC - NC[i][0] gives the x - coordinate of the ith global node, NC[i][1] gives the y
    // Here, 2 dofs make up one node and hence pairs of dofs will have the same coordinates
    for(int i = 0; i < total_dofs - 1 ; i = i + dim){
        NC[i][0] = x;
        NC[i+1][0]  = x; 
        x += incr_x;
        NC[i][1] = y;
        NC[i+1][1] = y;
        // If we have reached the x limit, reset x  to 0 and increment y
        if(abs(NC[i][0] - L) < 0.0000001){
            x = 0;
            y += incr_y;
        }

    }
    no_of_nodes_per_element = 4;
    // Since each node has more than 1 dof, the dofs per element will be the no of nodes per element * dimensions
    dofs_per_ele = no_of_nodes_per_element * dim;
    nel = nelx_ * nely_;
    EC.resize(nel);
    // Since EC is a vector of a vector , we need to initialize each row of EC
    for(int i = 0; i < nel; i++){
        // Over here we have to use dofs_per_ele as these will be the number of columns in EC
        EC[i] = std::vector<int>(dofs_per_ele);
    }
    nnx = nelx_ + 1; // Number of nodes along x
    nny = nely_ + 1; // Number of nodes along y
    unsigned int dofs_x = nnx * 2; // Number of dofs along x
    unsigned int dofs_y = nny * 2; // Number of dofs along y
    int inc_x = 0; // Tells how many increments we have had in x
    int inc_y = 0; // Tells how many increments we have had in y
    unsigned int n_count = 0;
    
    // Construct EC - EC[i][j] gives the global node number for local node j in element i
    for(int i = 0; i < nel;i++){
//            If we have reached last node on x, we increment y and reset our x coutner
        if(inc_x == dofs_x - 2){
            inc_y+=2;
            inc_x = 0;
        }
//      Storing clockwise on each element - pattern was hand derived using some examples for lesser number of elements
//		Node numbers increase left to right. Inc_y takes us to the nodes of the element 1 level down
        EC[i][0] = n_count + inc_y;
        EC[i][1] = EC[i][0] + 1;
        EC[i][2] = EC[i][1] + 1;
        EC[i][3] = EC[i][2] + 1;
        EC[i][4] = dofs_x + n_count + inc_y + 2; //Clock wise direction
        EC[i][5] = EC[i][4] + 1;
        EC[i][6] = dofs_x + n_count + inc_y;
        EC[i][7] = EC[i][6] + 1;
        inc_x += 2;
        n_count +=2;
    }
    // Set up quadrature data - Cant change number of quad points for now - Can include functionality with simple if-else if needed
    quad_rule = no_quad_points;
    quad_points.resize(quad_rule); //Resize quadpoints to appropriate size
    for(int i = 0; i < quad_rule; i++){
        quad_points[i] = std::vector<double>(dim);
    }
    quad_weights.resize(quad_rule); //Resize quadweights to appropriate size

    quad_points[0][0] = -sqrt(3./5.); // xi
    quad_points[0][1] = -sqrt(3./5.); // eta
    quad_points[1][0] = 0;
    quad_points[1][1] = 0;
    quad_points[2][0] = sqrt(3./5.);
    quad_points[2][1] = sqrt(3./5.);

    quad_weights[0] = 5./9.;
    quad_weights[1] = 8./9.;
    quad_weights[2] = 5./9.;

}

void FE::define_boundary_condition(double force, double g){
    std::cout<<std::endl<<"Defining boundary condition"<<std::endl;
//    Initialize the Dirichlet and Neumann boundary condtions
    g1 = g;
    f1 = force;
//    At each dof which is a boundary, we will put the value of g1, else we will put 0
    boundary_values.resize(total_dofs);
//    This function defines the boundary condition for a cantilivered beam i.e. all dofs at x = 0 have 0 displacement
    for(unsigned int dof_no = 0; dof_no < total_dofs ; dof_no++){
        if(NC[dof_no][0] == 0){
            boundary_values[dof_no] = g1;
            boundary_nodes.push_back(dof_no);
        }
    }
    // We define the F matrix fully here itself as we have no body force and just a force on the bottom right node acting downwards - Note, the way NC is set up, downwards is +ve Y axis and east is +ve x axis
    for(unsigned int dof_no = 0; dof_no < total_dofs ; dof_no++){
        if((abs(NC[dof_no][0] - L) < 0.00001) && (abs(NC[dof_no][1] - B) < 0.00001)){
            F[dof_no] = 0; // There are 2 dofs that satisfy this constraint - the first one is the x dof so this will be 0
            F[dof_no + 1] = f1; // The second dof is the y dof and this sbould have a force
            break; // After we have found this dof, we break otherwise we will be setting a force for someother dof
        }
    }
}


void FE::init_data_structs(){
    std::cout<<"Initializing data structures"<<std::endl;
    K.resize(total_dofs,total_dofs); //Resize K
    K.setZero(total_dofs,total_dofs); // Initialize K to 0
    F.resize(total_dofs); //Resize F
    F.setZero(total_dofs); // Setting F to zero here itself since we know the size
    U.resize(total_dofs); //Resive d
}

// Function for calculating the value of C - elasticity tensor

double FE::C(unsigned int i, unsigned int j, unsigned int k, unsigned int l){
    double lambda = (E * nu)/((1. + nu) * (1. - 2.*nu));
    double mu = E/(2. *(1. + nu));
    return lambda * (i==j) * (k==l) + mu * ((i==k)*(j==l) + (i==l)* (j==k));
}



void FE::fe_impl(Eigen::MatrixXd x){
    std::cout<<"Starting FE implementation"<<std::endl;
    // Initializing Flocal and Klocal
    std::vector<double> Flocal(dofs_per_ele);
    std::vector<std::vector<double> > Klocal(dofs_per_ele);
    for(int res = 0; res < dofs_per_ele; res++){
        Klocal[res] = std::vector<double>(dofs_per_ele);
    }
    
    // Initialize Jacobian matrix
    Eigen::MatrixXd Jac;
    Jac.resize(dim,dim);
    // Initialize Inverse Jacobian matrix
    Eigen::MatrixXd invJ;
    invJ.resize(dim,dim);
    // declare determinate of J
    double detJ;
    
    for(int ele = 0; ele < nel ; ele++){
        detJ = 0; // Set determinent of J to 0 - In this case it does not vary with each element, however for the generatl case, it should be within the element loop
        //      Initialize all elements in Klocal to 0
        std::fill(Klocal.begin(), Klocal.end(), std::vector<double>(dofs_per_ele, 0.));
        // Set all elements of Flocal to 0 - we will change the boundary elements to the force towards the end of this function
        std::fill(Flocal.begin(), Flocal.end() , 0.);
        
        for(unsigned int q1 = 0; q1 < quad_rule ; q1++){
            for(unsigned int q2 = 0; q2 < quad_rule ; q2++){
                Jac.setZero(dim,dim); // Reset J to zero for all quad points
                // Looping through the dimensions
                for(unsigned int i = 0; i < dim; i++){
                    for(unsigned int j = 0; j < dim; j++){
                        // Looping through the nodes of an element
                        for(unsigned int A = 0; A < no_of_nodes_per_element; A ++){
                            // Over here dim*A is used because EC has dim dofs per node. Each of these dofs have the same coordinate, so we can pick either one while calculating the jacobian. Over here, we use all the even dofs
                            Jac(i,j) += NC[EC[ele][dim*A]][i] * basis_gradient(A, quad_points[q1][i], quad_points[q2][j])[j];
                        }
                    }
                }
                detJ = Jac.determinant();
                invJ = Jac.inverse();
                // Now we go ahead and fill in the Klocal array
                for(unsigned int A = 0; A < no_of_nodes_per_element; A++){
                    // Capital I and K denote the physical coordinates
                    for(unsigned int I = 0; I < dim; I++){
                        for(unsigned int B=0 ; B < dim; B++){
                            for(unsigned int K = 0; K < dim; K++){
                                for(unsigned int J = 0; J < dim; J++){
                                    for(unsigned int L = 0; L < dim; L++){
                                        // Looping over the parametric coordinates - I think we only need to loop over j and k since only those indicies are used - Not sure though
                                        for(unsigned int j = 0; j < dim; j++){
                                            for(unsigned int l = 0; l < dim; l++){
                                                // Added i and k since we maybe do need it - Need to figure out how to reduce these number of loops - Will be too slow
                                                for(unsigned int i = 0 ; i < dim ; i++){
                                                    for(unsigned int k = 0; k < dim ; k ++){
                                                        Klocal[dim*A + I][dim*B + K] += (basis_gradient(A, quad_points[q1][0], quad_points[q2][1])[j] * invJ(j,J)) * C(I,J,K,L) * (basis_gradient(A, quad_points[q1][0], quad_points[q2][1])[l] * invJ(l,L)) * detJ * quad_weights[q1] * quad_weights[q2];
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        
                    }
                }
            }
        }
        // Now we assemble the Klocal into the K matrix which is the global matrix
        for(unsigned int I = 0; I < dofs_per_ele ; I++){
            for(unsigned int J = 0; J < dofs_per_ele ; J++){
                K(EC[ele][I],EC[ele][J]) += Klocal[I][J];
            }
        }
    }
    std::cout << "The determinant of K is " << K.determinant() << std::endl;
    
}
