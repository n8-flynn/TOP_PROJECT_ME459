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
    no_of_nodes = (nelx_ + 1) * (nely_ + 1);/**< Defines the total number of nodes in the domain*/
	std::cout<<"Total no. of nodes is "<<no_of_nodes<<std::endl;
    // Each node has 2 degrees of freedom as the number of dimensions is 2
    dim = 2;/**< Defines the dimension of the problem. The code only currently works for a dimension of 2.*/
    total_dofs = no_of_nodes * dim;/**< Since each node has 2 degrees of freedom (x and y), this variable defines the total number of degrees of freedom in the domain given by (no_of_nodes * dim)*/

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
    no_of_nodes_per_element = 4;/**<Number of nodes per element -  Since we are using Bi-linear Lagrange Basis functions, the number of nodes per element is always 4.*/
    // Agian , since each node has more than 1 dof, the dofs per element will be the no of nodes per element * dimensions
    dofs_per_ele = no_of_nodes_per_element * dim;/**< The degrees of freedoms per element -  Again the (number of nodes per element * dimension)*/
    nel = nelx_ * nely_; /**< Number of Elements - Just the product of number of elements along each direction i*/
    
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
    unsigned short int inc_x_ = 0; // Tells how many increments we have had in x
    unsigned short int inc_y_ = 0; // Tells how many increments we have had in y
    
    // Construct EC - EC[i][j] gives the global node number for local node j in element i
    for(unsigned short int i = 0; i < nel;i++){
//            If we have reached last node on x, we increment y and reset our x coutner
        if(inc_y_ == nny_ - 1){
            inc_x_+=1;
            inc_y_ = 0;
        }
//            Storing clockwise on each element - pattern was hand derived using some examples for lesser number of elements
//            Node numbers increase left to right. Inc_y takes us to the nodes of the element 1 level down
        EC_2[i][0] = i + inc_x_;
        EC_2[i][1] = i + 1 + inc_x_;
        EC_2[i][2] = nny_ + 1 + i + inc_x_;
        EC_2[i][3] = nny_ + i + inc_x_;
        inc_y_ += 1;
    }
    
    //unsigned short int nnx = nelx_ + 1; // Number of nodes along x
    unsigned short int nny = nely_ + 1; // Number of nodes along y
   // unsigned short int dofs_x = nnx * 2; // Number of dofs along x
    unsigned short int dofs_y = nny * 2; // Number of dofs along y
    unsigned short int inc_x = 0; // Tells how many increments we have had in x
    unsigned short int inc_y = 0; // Tells how many increments we have had in y
    unsigned short int n_count = 0;
    // Construct EC - EC[i][j] gives the global node number for local node j in element i
    for(unsigned short int i = 0; i < nel;i++){
   //            If we have reached last node on x, we increment y and reset our x coutner
           if(inc_y == dofs_y - 2){
               inc_x+=2;
               inc_y = 0;
           }
       //      Storing clockwise on each element - pattern was hand derived using some examples for lesser number of elements
       //        Node numbers increase left to right. Inc_y takes us to the nodes of the element 1 level down
               EC[i][0] = n_count + inc_x;
               EC[i][1] = EC[i][0] + 1;
               EC[i][2] = EC[i][1] + 1;
               EC[i][3] = EC[i][2] + 1;
               EC[i][4] = dofs_y + n_count + inc_x + 2; //Clock wise direction
               EC[i][5] = EC[i][4] + 1;
               EC[i][6] = dofs_y + n_count + inc_x;
               EC[i][7] = EC[i][6] + 1;
               inc_y += 2;
               n_count +=2;
    }
    // Set up quadrature data - Cant change number of quad points for now - Can include functionality with simple if-else if needed
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

void FE::define_boundary_condition(double force, double g,int wh){
    std::cout<<std::endl<<"Defining boundary condition"<<std::endl;
//    Initialize the Dirichlet and Neumann boundary condtions
    g1 = g;
    f1 = force;
//    At each dof which is a boundary, we will put the value of g1, else we will put 0
    boundary_values.resize(total_dofs);
//    This function defines the boundary condition for a cantilivered beam i.e. all dofs at x = 0 have 0 displacement
    for(unsigned short int dof_no = 0; dof_no < total_dofs ; dof_no++){
        if(NC[dof_no][0] == 0){
            boundary_values[dof_no] = g1;
            boundary_nodes.push_back(dof_no);
        }
        // We define the F matrix fully here itself as we have no body force and just a force on the bottom right node acting downwards - Note, the way NC is set up, downwards is +ve Y axis and east is +ve x axis
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

void FE::init_data_structs(){
    std::cout<<"Initializing data structures"<<std::endl;
    K.resize(total_dofs,total_dofs); //Resize K
    K.setZero();
    F.resize(total_dofs); //Resize F
    F.setZero(total_dofs); // Setting F to zero here itself since we know the size
    U.resize(total_dofs); //Resive d
}

// Function for calculating the value of C - elasticity tensor

inline double FE::C(uint8_t i, uint8_t j, uint8_t k, uint8_t l){
    double lambda = (E * nu)/((1. + nu) * (1. - 2.*nu));
    double mu = E/(2. *(1. + nu));
    return lambda * (i==j) * (k==l) + mu * ((i==k)*(j==l) + (i==l)* (j==k));
}

inline void FE::cal_jac(uint8_t q1, uint8_t q2){
    Eigen::MatrixXd Jac;
    Jac.resize(dim,dim);
    invJ.resize(dim,dim);
    for(uint8_t i = 0; i < dim; i++){
        for(uint8_t j = 0; j < dim; j++){
            Jac(i,j) = 0;
            // Looping through the nodes of an element
            for(unsigned short int A = 0; A < no_of_nodes_per_element; A++){
                // Over here dim*A is used because EC has dim dofs per node. Each of these dofs have the same coordinate, so we can pick either one while calculating the jacobian. Over here, we use all the even dofs
                Jac(i,j) += NC[EC[0][dim*A]][i] * basis_gradient(A, quad_points[q1][0], quad_points[q2][1])[j];
            }
        }
    }
    detJ = Jac.determinant();
    invJ = Jac.inverse();
//    saveData("invJ.csv", invJ);
}

void FE::cal_k_local(){
    std::cout<<"Determining Klocal"<<std::endl;
    // Initializing Klocal
    Klocal.resize(dofs_per_ele,dofs_per_ele);
    Klocal.setZero();
//    for(int res = 0; res < dofs_per_ele; res++){
//        Klocal[res] = std::vector<double>(dofs_per_ele);
//    }
    

//    std::fill(Klocal.begin(), Klocal.end(), std::vector<double>(dofs_per_ele, 0.));
    for(uint8_t q1 = 0; q1 < quad_rule ; q1++){
        for(uint8_t q2 = 0; q2 < quad_rule ; q2++){
            cal_jac(q1,q2);
            // Now we go ahead and fill in the Klocal array
            for(unsigned short int A = 0; A < no_of_nodes_per_element; A++){
                // Capital I and K denote the physical coordinates
                for(uint8_t I = 0; I < dim; I++){
                    for(unsigned short int B=0 ; B < no_of_nodes_per_element; B++){
                        for(uint8_t K = 0; K < dim; K++){
                            for(uint8_t J = 0; J < dim; J++){
                                for(uint8_t L = 0; L < dim; L++){
                                    // Looping over the parametric coordinates - I think we only need to loop over j and k since only those indicies are used - Not sure though
                                    for(uint8_t j = 0; j < dim; j++){
                                        for(uint8_t l = 0; l < dim; l++){
                                            // Added i and k since we maybe do need it - Need to figure out how to reduce these number of loops - Will be too slow
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
//    saveData("k_local.csv", Klocal);
}


void FE::assemble(Eigen::MatrixXd x,double penal){
    //std::cout<<"Assembling and applying dirichlet conditions"<<std::endl;
//    Klocal<<0.4945,0.1786,-0.3022,-0.0137,-0.2473,-0.1786,0.0549,0.0137,
//            0.1786,0.4945,0.0137,0.0549,-0.1786,-0.2473,-0.0137,-0.3022,
//            -0.3022,0.0137,0.4945,-0.1786,0.0549,-0.0137,-0.2473,0.1786,
//            -0.0137,0.0549,-0.1786,0.4945,0.0137,-0.3022,0.1786,-0.2473,
//            -0.2473,-0.1786,0.0549,0.0137,0.4945,0.1786,-0.3022,-0.0137,
//            -0.1786,-0.2473,-0.0137,-0.3022,0.1786,0.4945,0.0137,0.0549,
//            0.0549,-0.0137,-0.2473,0.1786,-0.3022,0.0137,0.4945,-0.1786,
//            0.0137,-0.3022,0.1786,-0.2473,-0.0137,0.0549,-0.1786,0.4945;
    unsigned short int ely = 0;
    unsigned short int elx = 0;
    double x_;
//    These are all dummy variables
    unsigned short int row1;
    unsigned short int col1;
    unsigned short int row2;
    unsigned short int col2;
    for(unsigned short int ele = 0; ele < nel ; ele++){
        if(ely == nely_){
            ely = 0;
            elx++;
        }
        x_ = x(ely,elx);
        ely++;
        // Now we assemble the Klocal into the K matrix which is the global matrix
        for(unsigned short int I = 0; I < dofs_per_ele ; I++){
            row1 = EC[ele][I];
            row2 = I;
            for(unsigned short int J = 0; J < dofs_per_ele ; J++){
                col1 = EC[ele][J];
                col2 = J;
                K(row1, col1) += pow(x_,penal) * Klocal(row2,col2);
                
            }
        }
    }
//    std::cout << "The determinant of K is " << K.determinant() << std::endl;
//    saveData("k_before.csv", K);
    // Now we apply the Dirichlet boundary conditons and modify K accordingly
   // std::cout<<"Applying Dirichlet BC's"<<std::endl;
    for(unsigned short int i : boundary_nodes){
        double g = boundary_values[i];
        // Loop to move the approprate column of K to the RHS - source - https://www.math.colostate.edu/~bangerth/videos.676.21.65.html
        for(unsigned short int row = 0; row < total_dofs; row++){
            // This condition is so that a dof which has already been set in F is not changed
            if(row == i){
                continue;
            }
            // All the other dofs in F are varied as we move the column of K to the RHS
            else{
                F[row] = F[row] - g * K(row,i);
            }
        }
        // Set all the diagonal elements to 1 and all the other elements in the row and column to 1
        K.row(i) *= 0;
        K.col(i) *= 0;
        K(i,i) = 1.;
        // Set the value in F at athe node
        F[i] = g;
    }
//    std::cout << "After applying BC the determinant of K is " << K.determinant() << std::endl;
//    saveData("k_after.csv", K);
//    std::cout<<K<<std::endl;
}
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
