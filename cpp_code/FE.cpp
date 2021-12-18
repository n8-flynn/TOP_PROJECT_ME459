// Created by Huzaifa Mustafa Unjhawala 

#include "FE.h"
#include <fstream>
#include <iostream>

/** \brief The FE Class constructor.
 * 
 * \param nelx Defines the number of elements along the X-direction.
 * \param nely Defines the number of elements along the Y-direction.
 * \param length Defines the length (distance along x). However, in the current version of the code the length always equals nelx.
 * \param breadth Defines the breadth (distance along y). However, in the current version of the code the length always equals nely.
 * \param youngs_mod Defines the Youngs Modulus (E) of the material.
 * \param pois_rat Defines the Poisson Ratio \f$(\nu)\f$ of the material.
 */

FE::FE(unsigned short int nelx, unsigned short int nely, double length, double breadth, double youngs_mod, double pois_rat){
	L = length;
	B = breadth;
	nelx_ = nelx;
	nely_ = nely;
	E = youngs_mod;
	nu = pois_rat;
}

// Basis function - Internal function needed for fe implementation
/**\brief The function to define the Bi-linear Lagrange Basis Functiions.
 *Calculates the value of the Bi-linear basis functions corresponding to the node "node"  at the point 'xi' and 'eta' in the parametric space.
 *The function is made in-line for better performance.
 * \param node This is the node to which the basis function belongs to.
 * \param xi The \f$(\xi)\f$ component of the point in the parametric space where we are evaluating the basis function of a given node.
 * \param eta The \f$(\eta)\f$ component of the point in the parametric space where we are evaluating the basis function of a given node.
 */

inline double FE::basis_function(unsigned short int node , double xi, double eta){
//Here I use a switch case to identify the node we and then apply the xi and eta value to evaluate the output. The default case will send out an error message. However, the program will not stop. Exception needs to be added.
    double output;
    switch(node) {
        case 0:
            output = 0.25 * (1-xi) * (1-eta);
            break;
        case 1:
            output = 0.25 * (1+xi) * (1-eta);
            break;
        case 2:
            output = 0.25 * (1+xi) * (1+eta);
            break;
        case 3:
            output = 0.25 * (1-xi) * (1+eta);
            break;
        default:
            std::cout<<"There is no "<<node<<" node, invalid input"<<std::endl;
            break;
    }
	return output;
}
// Basis function defined using the general formula obtained from

/** \brief The function to calculate the gradient of the Bi-linear Lagrange Basis Function.
 * Calculates the gradient of the Bi-linear basis function corresponding to the node "node" and at the point 'xi' and 'eta' with respect to both 'xi' and 'eta'.
 * Thus, this method returns a vector of length 2 whose 0th element is the gradient with respect to 'xi' and the 1st element is the gradient with respect to 'eta'
 * \param node This is the node to which the basis function whose gradient we are taking belongs to
 * \param xi The \f$(\xi)\f$ component of the point in the parametric space where we are evaluating the gradient of the basis function of a given node.
 * \param eta The \f$(\eta)\f$ component of the point in the parametric space where we are evaluating the gradient of the basis function of a given node.
 */
inline std::vector<double> FE::basis_gradient(unsigned short int node,double xi, double eta){
// Similar logic to basis_function.he default case will send out an error message. However, the program will not stop. Exception needs to be added.
    
    std::vector<double> bg(dim,0.0);
    switch(node) {
        case 0:
            bg[0] = -0.25 * (1 - eta);
            bg[1] = -0.25 * (1 - xi);
            break;
        case 1:
            bg[0] = 0.25 * (1 - eta);
            bg[1] = -0.25 * (1 + xi);
            break;
        case 2:
            bg[0] = 0.25 * (1 + eta);
            bg[1] = 0.25 * (1 + xi);
            break;
        case 3:
            bg[0] = -0.25 * (1 + eta);
            bg[1] = 0.25 * (1 - xi);
            break;
        default:
            std::cout<<"There is no "<<node<<" node, invalid input"<<std::endl;
            break;
    }

    return bg;
}

/** \brief The Mesh generator.
 * This method fills up the Nodal Coordinate matrix and the Elemental Connectivity matrix whihc define the mesh of the domain.
 * It also takes as input the number of quadrature points in order to define the location of the quadrature points and the quadrature weights.
 * Currently this method can only take 3 quadrature points, however, this can easily be explanded in the future.
 * \param no_quad_points Defines the number of Guassian Quadrature Points taken along each direction in the parametric space.
*/
void FE::mesh(uint8_t no_quad_points){
    
    // Each element is made of of 2 nodes. Thus along each direction we will have (number of elements + 1) nodes.
    no_of_nodes = (nelx_ + 1) * (nely_ + 1);// Defines the total number of nodes in the domain
	std::cout<<"Total no. of nodes is "<<no_of_nodes<<std::endl;
    // Each node has 2 degrees of freedom as the number of dimensions is 2
    dim = 2;// Defines the dimension of the problem. The code only currently works for a dimension of 2.
    total_dofs = no_of_nodes * dim;// Since each node has 2 degrees of freedom (x and y), this variable defines the total number of degrees of freedom in the domain given by (no_of_nodes * dim)

	// Nodal Coordinate array will give the coordinates of each dof not just each node
	NC.resize(total_dofs);
    // Since NC is a vector of a vector , we need to initialize each row of EC
    for(unsigned int i = 0; i < total_dofs; i++){
        NC[i] = std::vector<double>(dim,0.0);
    }



    // As I loop through all the dofs, incr_x helps me identify the which dof I am on along the x direction. Since all the nodes are equally spaced, the incr_x will be constant as well move along the x direction.
    
	double incr_x = L/(nelx_);
    
    // As I loop through all the dofs, incr_y helps me identify the which dof I am on along the y direction. Since all the nodes are equally spaced, the incr_y will be constant as well move along the y direction.The incr_y is negetive because we start from the top left corner which has the coordinates (0,B) and then move downwards. Hence, incr_y will always be negetive.
    
    double incr_y = -B/(nely_);
    
    // First node is at (0,B)
    double x = 0.0;
    double y = B;
    
    // I now Construct NC - NC[i][0] gives the x - coordinate of the ith global node, NC[i][1] gives the y
    // Here, 2 dofs make up one node and hence pairs of dofs will have the same coordinates
    // Even thoough 'dimension' number of dofs (rows in NC) will have the same values, it is important to define NC in this way as it makes the loop to calculate the Elemental Stiffnes Matrix (Klocal) much easier.
    
    // We increment by dim as we set both dofs of a paticular node to their respective coordinate. Makes the code a bit faster with loop unrolling.
    for(unsigned short int i = 0; i < total_dofs - 1 ; i = i + dim){
        NC[i][0] = x;
        NC[i+1][0]  = x; 
        
        NC[i][1] = y;
        NC[i+1][1] = y;
        // We move one step down towards (0,0) along the y axis. Thus we only need to change the y coordiante until we hit the last element along the y direction.
        
        y += incr_y;
        
        // Check if we y coordinate of our dof is at 0. Basically tells if we have reached all the way down and need to shift to the right and go back to the top of the domain.
        if(abs(NC[i][1]) < 0.0000001){
            // Increment x and set y again to the top - breadth.
            x += incr_x;
            y = B;
        }

    }
    // Note - Have to decide whether to remove EC_2 or keep it
    
    no_of_nodes_per_element = 4;//<Number of nodes per element -  Since we are using Bi-linear Lagrange Basis functions, the number of nodes per element is always 4.
    // Agian , since each node has more than 1 dof, the dofs per element will be the no of nodes per element * dimensions
    dofs_per_ele = no_of_nodes_per_element * dim;//< The degrees of freedoms per element -  Again the (number of nodes per element * dimension)
    nel = nelx_ * nely_; //< Number of Elements - Just the product of number of elements along each direction i
    
    EC_2.resize(nel);
    EC.resize(nel);
    // Since EC is a vector of a vector , we need to initialize each row of EC
    for(unsigned short int i = 0; i < nel; i++){
        // Over here we have to use dofs_per_ele as these will be the number of columns in EC
        EC[i] = std::vector<unsigned short int>(dofs_per_ele);
        EC_2[i] = std::vector<unsigned short int>(no_of_nodes_per_element);
    }
    
    
    nnx_ = nelx_ + 1; // Number of nodes along x
    nny_ = nely_ + 1; // Number of nodes along y
    unsigned short int inc_x_ = 0; //Again, like in NC - tells how many increments we have had in x
    unsigned short int inc_y_ = 0; // Again like in NC - tells how many increments we have had in y
    
    // Construct EC_2 - EC_2[i][j] gives the global node number for local node j in element i. Stress on node number as that is the difference between EC and EC_2. EC_2 gives the global node number wheras EC gives the gloabl degree of freedom.
    for(unsigned short int i = 0; i < nel;i++){
        // Similar to NC, if we reach the last element along the Y direction for a paticular X, we increment X and set the uncr_y to 0. Basically moving to one column on the right
        if(inc_y_ == nny_ - 1){
            inc_x_+=1;
            inc_y_ = 0;
        }
        // The local node numbering is in the anti clockwise direction as can be seen below
        // 1---------4
        // -----------
        // -----------
        // -----------
        // 2---------3
        
        // So,for the first element this will correspond to the following global node numbers (nny is the number of nodes along y)
        
        //1-----------nny+1
        //----------------
        //----------------
        //----------------
        //2----------nny+2
        
        // Going off the above diagrams we get the following pattern. Loop unrolling done for simplicity of code and better performance.
        EC_2[i][0] = i + inc_x_;
        EC_2[i][1] = i + 1 + inc_x_;
        EC_2[i][2] = nny_ + 1 + i + inc_x_;
        EC_2[i][3] = nny_ + i + inc_x_;
        inc_y_ += 1;
    }
    
    unsigned short int dofs_y = nny_ * 2; // Number of dofs along y - needed as now we are constructing EC which gives the global dof number for a local dof number of a paticular element.
    // I set the our increments back to 0
    inc_x_ = 0;
    inc_y_ = 0;
    // I introduce a new variable n_count in order to count the total number of dofs we have completed accounting for. This was introduced to help with the algorithim as I could not think of any other way to do it.
    unsigned short int n_count = 0;
    // Construct EC - EC[i][j] gives the global dof number for local dof j in element i
    
    // Again, like above we go in an anticlockwise direction for our local dof numbering
    //1,2-------7,8
    //-----------
    //-----------
    //-----------
    //3,4------5,6
    
    // This corresponds to a global DOF of
    //1,2-------dofs_y+1,dofs_y+2
    //----------
    //----------
    //----------
    //3,4-------dofs_y+3,dofs_y+4
    
    
    for(unsigned short int i = 0; i < nel;i++){
            //
           if(inc_y_ == dofs_y - 2){
               inc_x_+=2;
               inc_y_ = 0;
           }
                // Again going off the above logic, this is the best way I could think off doing it. Again, there is sone loop unrolling and uses some elements that are already present in the cache (temporal and spatial locality) for better performance.
               EC[i][0] = n_count + inc_x_;
               EC[i][1] = EC[i][0] + 1;
               EC[i][2] = EC[i][1] + 1;
               EC[i][3] = EC[i][2] + 1;
               EC[i][4] = dofs_y + n_count + inc_x_ + 2;
               EC[i][5] = EC[i][4] + 1;
               EC[i][6] = dofs_y + n_count + inc_x_;
               EC[i][7] = EC[i][6] + 1;
               inc_y_ += 2;
               n_count +=2;
    }

    // The below code defines the Guassian Quadrature and weights.
    
    quad_rule = no_quad_points;
    quad_points.resize(quad_rule); //Resize quadpoints to appropriate size
    for(uint8_t i = 0; i < quad_rule; i++){
        quad_points[i] = std::vector<double>(dim);
    }
    quad_weights.resize(quad_rule); //Resize quadweights to appropriate size

    quad_points[0][0] = sqrt(3./5.); // xi
    quad_points[0][1] = sqrt(3./5.); // eta
    quad_points[1][0] = 0;
    quad_points[1][1] = 0;
    quad_points[2][0] = -sqrt(3./5.);
    quad_points[2][1] = -sqrt(3./5.);

    quad_weights[0] = 5./9.;
    quad_weights[1] = 8./9.;
    quad_weights[2] = 5./9.;
}
/** \brief Defines the boundary conditions on the problem.
 * Given the force, the location of the force and the dirichlet displacement, this function assembles the boundary values and the boudary nodes vectors. The F vector which makes up the RHS of the system we eventually solve is also modified with the force value added to the appropriate row.
 * \param force The value of the external force that we apply on our domain
 * \param g The value of the dirichlet boundary condition (displacement)
 * \param wh This tells us where the force is to be applied. A value of 0 applies the force at (L,0). A value of 1 applies the force at (L/2,B/2). A value of 2 applied the force at (L,B).
 */
void FE::define_boundary_condition(double force, double g,int wh){
    std::cout<<std::endl<<"Defining boundary condition"<<std::endl;
    // Assign the required inputs to class members
    g1 = g;
    f1 = force;
//    At each dof which is a boundary, I will put the value of g1, else we will put 0
    boundary_values.resize(total_dofs);
//    This function defines the boundary condition for a cantilivered beam i.e. all dofs at x = 0 have 0 displacement
    for(unsigned short int dof_no = 0; dof_no < total_dofs ; dof_no++){
        if(NC[dof_no][0] == 0){
            boundary_values[dof_no] = g1;
            boundary_nodes.push_back(dof_no);
        }
        // We assign the force to the F vector here - Note, the way NC is set up, downwards is +ve Y axis and east is +ve x axis. Another note - There is no body force and all the dirichlet boundary conditions are 0 and hence this will be all there is to the F vector.
        // Based on wh, we apply the force to the appropriate node and DOF (Y axis in our case as we are applying a downward force).
        if(wh == 0){
            if((abs(NC[dof_no][0] - L) < 0.00001) && (abs(NC[dof_no][1]) < 0.00001)){
                F[dof_no] = 0; // There are 2 dofs that satisfy this constraint - the first one is the x dof so this will be 0
                F[dof_no + 1] = -f1; // The second dof is the y dof and this sbould have a force
                break; // After we have found this dof, we break otherwise we will be setting a force for someother dof
            }
        }
        else if(wh == 1){
            if((abs(NC[dof_no][0] - L) < 0.00001) && (abs(NC[dof_no][1] - B/2) < 0.00001)){
                F[dof_no] = 0; // There are 2 dofs that satisfy this constraint - the first one is the x dof so this will be 0
                F[dof_no + 1] = -f1; // The second dof is the y dof and this sbould have a force
                break; // After we have found this dof, we break otherwise we will be setting a force for someother dof
            }
        }
        else{
            if((abs(NC[dof_no][0] - L) < 0.00001) && (abs(NC[dof_no][1] - B) < 0.00001)){
                F[dof_no] = 0; // There are 2 dofs that satisfy this constraint - the first one is the x dof so this will be 0
                F[dof_no + 1] = -f1; // The second dof is the y dof and this sbould have a force
                break; // After we have found this dof, we break otherwise we will be setting a force for someother dof
            }
        }
    }

}
/** \brief A method to save Eigen Matrices and Vectors as CSV files for debugging purposes.
 * \param fileName The name of the csv file to be saved
 * \param matrix The Eigen Matrix that is to be stored as a csv file
 */
void FE::saveData(std::string fileName, Eigen::MatrixXd  matrix)
{
    //https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
 
    std::ofstream file(fileName);
    if (file.is_open())
    {
        file << matrix.format(CSVFormat);
        file.close();
    }
}

/** \brief Initializes all the datastructres that are needed to solve the linear system of equations and deterimine U
 * The K matrix (Global Stiffness Matrix) is resized to size (total dofs, total dofs) and all elements are set to zero.
 * The F Vector is resized to size (total dofs) and all elements are set to zero.
 * The solution vector U (displacement vector) is resized to size (total dofs).
 */
void FE::init_data_structs(){
    std::cout<<"Initializing data structures"<<std::endl;
    K.resize(total_dofs,total_dofs); //Resize K
    K.setZero();
    F.resize(total_dofs); //Resize F
    F.setZero(total_dofs); // Setting F to zero here itself since we know the size
    U.resize(total_dofs); //Resive d
}

/** \brief This is an inline function used to compute the value of the elasticity tensor at i,j,k and l.
 * I am using a standard linear isotropic material and thus \f$(\lambda = \frac{\nu E}{(1 + \nu)(1-2\nu)})\f$ and \f$(\mu = \frac{E}{2(1+\nu)})\f$.
 * C is thus \f$(\lambda\delta_{IJ}\delta_{KL} + \mu(\delta_{IK}\delta_{JL} + \delta_{IL}\delta_{JK}))\f$. Where \f$(\delta)\f$ is the kronecker detla.
 * \param i Pyhsical dimension that goes from 0 to dimension - 1.
 * \param j Pyhsical dimension that goes from 0 to dimension - 1.
 * \param k Pyhsical dimension that goes from 0 to dimension - 1.
 * \param l Pyhsical dimension that goes from 0 to dimension - 1.
 */
// Function for calculating the value of C - elasticity tensor

inline double FE::C(uint8_t i, uint8_t j, uint8_t k, uint8_t l){
    double lambda = (E * nu)/((1. + nu) * (1. - 2.*nu));
    double mu = E/(2. *(1. + nu));
    return lambda * (i==j) * (k==l) + mu * ((i==k)*(j==l) + (i==l)* (j==k));
}


/** \brief This function is used to evaluate the Jacobian matrix, its inverse and its determinant at a paticular guassian quadrature point.
 * The Jacobian is calculated in a speerate inline function to prevent clutter in the already cluttered elemental K code.
 * \param q1 Identifies the number of the quadrature point along the 'xi' direction
 * \param q2 Identifies the number of the quadrature point along the 'eta' direction
 */
inline void FE::cal_jac(uint8_t q1, uint8_t q2){
    Eigen::MatrixXd Jac; // Define Jacobian locally as we only really need inverse of jacobian
    Jac.resize(dim,dim); // Resize jacobian to appropirate size
    invJ.resize(dim,dim); // Resize our class member Inverse of Jacobian to appropriate size
    for(uint8_t i = 0; i < dim; i++){
        for(uint8_t j = 0; j < dim; j++){
            Jac(i,j) = 0;
            // Looping through the nodes of an element
            for(unsigned short int A = 0; A < no_of_nodes_per_element; A++){
                // Over here dim*A is used because EC has dim dofs per node. Each of these dofs have the same coordinate, so we can pick either one while calculating the jacobian. Over here, we use all the even dofs
                // quad_points[q1][0] gives the 'xi' component of quad point q1 wheras using quad_points[q2][1] gives us the 'eta component.
                Jac(i,j) += NC[EC[0][dim*A]][i] * basis_gradient(A, quad_points[q1][0], quad_points[q2][1])[j];
            }
        }
    }
    // Standard methods on Eigen Matrices to get the determinant and the inverse
    detJ = Jac.determinant();
    invJ = Jac.inverse();
}


/**\brief This method is used to fill up the elemental K matrix.
 * Since our domain is made up of all the same type of elements with equal sizing, the cal_k_local method is only called once to fing the K elemental for the 1st element since all the elements will have the same k local.
 */


void FE::cal_k_local(){
    std::cout<<"Determining Klocal"<<std::endl;
    // Resizing Klocal and setting all its elements to 0.
    Klocal.resize(dofs_per_ele,dofs_per_ele);
    Klocal.setZero();
    // First we loop over the quadrature points. Since we have quadrature points in both directions, we loop over the quad points 2 times to get the 9 totoal combinations of quadrature points.
    for(uint8_t q1 = 0; q1 < quad_rule ; q1++){
        for(uint8_t q2 = 0; q2 < quad_rule ; q2++){
            // For a patricular quad point, we calculate the jacobian. This will fill up the inverse Jacobian matrix and calculate the determinant of the jacobian.
            cal_jac(q1,q2);
            // Then we loop over the nodes in order to recreate the weak form of the problem - we first loopover A for spacial locality
            for(unsigned short int A = 0; A < no_of_nodes_per_element; A++){
                // Capital I and K denote the physical coordinates - We first loop over I for spcial locality
                for(uint8_t I = 0; I < dim; I++){
                    for(unsigned short int B=0 ; B < no_of_nodes_per_element; B++){
                        for(uint8_t K = 0; K < dim; K++){
                            for(uint8_t J = 0; J < dim; J++){
                                for(uint8_t L = 0; L < dim; L++){
                                    // Looping over the parametric coordinates j and l.
                                    for(uint8_t j = 0; j < dim; j++){
                                        for(uint8_t l = 0; l < dim; l++){
                                            // Since we are only looping over the nodes per element, we need to use dim*A + I in order to fill all the degrees of freedom per element. Each evaluation requires alot of function calls, however, all of these functions are defined as inline functions and is thus hopefully faster.
                                            Klocal(dim*A + I,dim*B + K) += (basis_gradient(A, quad_points[q1][0], quad_points[q2][1])[j] * invJ(j,J)) * C(I,J,K,L) * (basis_gradient(B, quad_points[q1][0], quad_points[q2][1])[l] * invJ(l,L)) * detJ * quad_weights[q1] * quad_weights[q2];
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
/** \brief This method takes the K local found by cal_k_local and assembles the K Global Matrix.
 * In addition it also modifies the K global to account for the Dirichlet Boundary conditions. While assembling the K global matrix, the K local value at the appropriate local node is multiplied by the relative density 'x' of that paticular element raised to a power 'penal' which is the panalization power. In mathematical terms -
 * \f$(K(global dof A, global dof B) = x(global element number along y, global element number along x)^{penal} klocal(local dof A, local dof B))\f$.
 * Note - This is the function that is evaluated at every optimization loop and it is thus the most cruicial function in terms of performance.
 * \param x This is the relative densities Eigen Matrix that is used for the topology optimization.
 * \param penal The penalization power. The penalization power is used to refine the solution to solid and void regions to aid manufacturibility.As we can see from the formula, it will make the values closer to 0 go faster towards 0 and the values closer to 1 go faster towards 1.
 */

void FE::assemble(Eigen::MatrixXd x,double penal){
    // We need ely and elx as x is defined as a 2 d matrix where the relative density at (ely,elx) is the relative density at ely th row and the elx th column.
    unsigned short int ely = 0;
    unsigned short int elx = 0;
    double x_;
//    These are all dummy variables to have temporal locaility - Unsure if the compiler will take care of this
    unsigned short int row1; // Corresponds to the global node's row in K
    unsigned short int col1; // Corresponds to the global node's column in K
    unsigned short int row2; // Corresponds to the local node's row in klocal
    unsigned short int col2; // Corresponds to the local nodes column in klocal
    for(unsigned short int ele = 0; ele < nel ; ele++){
        if(ely == nely_){
            ely = 0;
            elx++;
        }
        // Evaluate it and store it here itself for temporal locality
        x_ = x(ely,elx);
        ely++;
        // Now we assemble the Klocal into the K matrix which is the global matrix
        // For spatial locality we loop over I first
        for(unsigned short int I = 0; I < dofs_per_ele ; I++){
            // We store row1 and row 2 as its not going to change in the next loop and it is costly evaluating it in each loop
            row1 = EC[ele][I];
            row2 = I;
            for(unsigned short int J = 0; J < dofs_per_ele ; J++){
                col1 = EC[ele][J];
                col2 = J;
                K(row1, col1) += pow(x_,penal) * Klocal(row2,col2);
                
            }
        }
    }
    // Now we apply the Dirichlet boundary conditons and modify K accordingly
    
    // First we loop over all the boundary nodes that were identified in the define_boundary_cond() function
    for(unsigned short int i : boundary_nodes){
        // Get the value g (displacement) at that boundary node - In our case, g is always 0
        double g = boundary_values[i];
        // Loop to move the approprate column of K to the RHS - The source of this algorithm is from - https://www.math.colostate.edu/~bangerth/videos.676.21.65.html
        for(unsigned short int row = 0; row < total_dofs; row++){
            // This condition is so that a dof which has already been set in F is not changed
            if(row == i){
                continue;
            }
            // All the other dofs in F are varied as we move the column of K to the RHS - In our case, F is not modified at all since g is zero. This loop can be removed if we seek more performance. However, for non zero g, this loop is necessary.
            // This basically moves the known varibales from the LHS to the RHS
            else{
                F[row] = F[row] - g * K(row,i);
            }
        }
        // Set all the diagonal elements to 1 and all the other elements in the row and column to 0
        K.row(i) *= 0;
        K.col(i) *= 0;
        K(i,i) = 1.;
        // Set the value in F at at the node - In our case it will always set to 0.
        F[i] = g;
    }
}


/** \brief Solves the linear system of equations KU = F using a Sparse LDLT Cholesky Decomoposition method provided by Eigen
 * The K and F global matrices that were defined and filled as Dense matrices, they are converted into sparse matrix using the sparseView() method available in Eigen. This is done because these matrices are in fact sparse and using sparse matrix over dense matrices dramatically speeds up the solving of the system of linear equations.
 * This method returns the solution vector U to the toplogy function where it is used to find the compliance 'c'.
 */
Eigen::VectorXd FE::solve(){
    //std::cout<<"."<<std::endl;
    Eigen::SparseMatrix<double> K_ = K.sparseView();
    Eigen::SparseVector<double> F_ = F.sparseView();
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    solver.compute(K_);
    U = solver.solve(F_);
//    std::cout<<U<<std::endl;
    return U;
}

/** \brief Used to write the solution U to a .vtu file for visulisation and debugging.s*/
void FE::fem_to_vtk(){
        
    // Write to file all the stuff that is needed for plotting
    std::cout<<"Writing to vtu file for steady state"<<std::endl;
    std::ofstream out_file;
    std::stringstream ss;
    ss<<"trial.vtu";
    std::string f_name;
    f_name+= ss.str();
//        out_file.open("res/output_" + std::to_string(scheme) + std::to_string(t_step) +  ".vtu");
    out_file.open(f_name);
    if(out_file.fail()){
        std::cout<<"File did not open"<<std::endl;
    }
//   Writing headers
    out_file<<"<?xml version=\"1.0\"?>"<<std::endl;
    out_file<<"<VTKFile type=\"UnstructuredGrid\"  version=\"0.1\"  >"<<std::endl;
    out_file<<"<UnstructuredGrid>"<<std::endl;
//    Writing nodal coordinates array
    out_file<<"<Piece  NumberOfPoints=\""<<no_of_nodes<<"\" NumberOfCells=\""<<nel<<"\">"<<std::endl;
    out_file<<"<Points>"<<std::endl;
    out_file<<"<DataArray type=\"Float32\" NumberOfComponents=\""<<dim+1<<"\" format=\"ascii\">"<<std::endl;
    float z = 0.;
    for(int node = 0; node < total_dofs; node = node+2){
        out_file<<(float) NC[node][0]<<" "<<(float) NC[node][1]<<" "<<z<<std::endl;
    }
    out_file<<"</DataArray>"<<std::endl;
    out_file<<"</Points>"<<std::endl;
    
//    Writing Element connectivity
    out_file<<"<Cells>"<<std::endl;
    out_file<<"<DataArray  type=\"UInt32\"  Name=\"connectivity\"  format=\"ascii\">"<<std::endl;
    for(int ele = 0; ele < nel; ele++){
        out_file<<EC_2[ele][0]<<" "<<EC_2[ele][1]<<" "<<EC_2[ele][2]<<" "<<EC_2[ele][3]<<std::endl;
    }

    out_file<<"</DataArray>"<<std::endl;
//    Writing element offsets vector(required in VTK format)
    unsigned int offsets = 0;
    out_file<<"<DataArray  type=\"UInt32\"  Name=\"offsets\"  format=\"ascii\">"<<std::endl;

    for(int ele = 0; ele < nel; ele++){
        offsets=offsets+4;
        out_file<<offsets<<std::endl;
    }


//    Writing element types vector(required in VTK format)
    out_file<<"</DataArray>"<<std::endl;
    out_file<<"<DataArray  type=\"UInt8\"  Name=\"types\"  format=\"ascii\">"<<std::endl;
    unsigned int ty = 9;
    for(int ele = 0; ele < nel; ele++){
        out_file<<ty<<std::endl;
    }
    out_file<<"</DataArray>"<<std::endl;
    out_file<<"</Cells>"<<std::endl;
//    Writing field values (if any) array
    
    out_file<<"<PointData  Scalars=\"u\">"<<std::endl;
    out_file<<"<DataArray  type=\"Float32\"  Name=\"Displacement\" NumberOfComponents=\"3\" format=\"ascii\">"<<std::endl;
    float z1 = 0;
    for(int node = 0; node < total_dofs; node=node+2){
        out_file<<(float)U(node)<<" "<<(float)U(node+1)<<" "<<z1<<std::endl;
//        out_file<<(float)d_n(node)<<std::endl;
    }
    out_file<<"</DataArray>"<<std::endl;
    out_file<<"</PointData> "<<std::endl;
    out_file<<"</Piece> "<<std::endl;
    out_file<<"</UnstructuredGrid> "<<std::endl;
    out_file<<"</VTKFile> "<<std::endl;
    out_file.close();
        
}
