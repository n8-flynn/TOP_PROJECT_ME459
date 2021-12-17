// Created by Huzaifa Mustafa Unjhawala

/** \file The FE Class header file
 */
#ifndef FE_H
#define FE_H

// All the required libraries
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
#include "Eigen"
//#include <Eigen/Dense>



typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;

/*! \brief The FE class computes the Stiffness Matrix and sovles for the displacement.
 *
 * The FE class takes as input the domain dimensions, mesh properties and the material properties.
 * It then computes the Global Stiffness Matrix (K), Global Force Vector (F) and computes the
 * Global Displacement (U) by solving the linear system KU = F.  Linear Lagrange Basis functions
 * are used to approximate U over the elements. A nine point Guass Quadrature is used to numerically
 * integrate the resulting equations. For more details, pleaase view documentation for class members and methods.
 */

class FE{
    public:
        FE(unsigned short int nelx,unsigned short int nely,double length, double breadth,double youngs_mod, double pois_rat);
        double basis_function(unsigned short int node, double xi,double eta);
        std::vector<double> basis_gradient(unsigned short int node, double xi,double eta);
        void mesh(uint8_t no_quad_points);
        void define_boundary_condition(double force, double g,int wh); // Function to fill the boundary_values (stores the values at the boundaries) and the boundary_nodes (stores the global node number of the nodes on the boundary)
        double C(uint8_t i, uint8_t j, uint8_t k, uint8_t l); // Used to get the elasticity tensor C
        void init_data_structs(); //To resize all the global matrices based on the mesh - Internal function
        void cal_jac(uint8_t q1, uint8_t q2);
        void cal_k_local(); // Calculates the K local for one element - As all elements are the same, can use the same klocal
        void assemble(Eigen::MatrixXd x,double penal); //Uses the klocal to assemble K global using the volume fractions
        Eigen::VectorXd solve(); // Solves and returns U which is then used in the toplogy code
        void fem_to_vtk();
        void saveData(std::string fileName, Eigen::MatrixXd  matrix);


        // Class datastructures
        double L,B,g1,f1,E,nu,lambda,mu,penal_,detJ; //Standard constants - L - Length, B - breadth, g1 - Dirichlet conditon, E - Youngs Modulus, nu - Poissons ration, lambda and mu
        // are the lame's parameters
        
        unsigned short int nelx_;/**< Number of elements along the x direction*/
        unsigned short int nely_;/**< Number of elements along the y direction*/
        unsigned short int nel;/**< Total Number of elements*/
        unsigned short int nnx_;/**< Number of nodes along the x direction*/
        unsigned short int nny_;/**< Number of nodes along the y direction*/
        unsigned short int no_of_nodes;/**< Total Number of nodes*/
        unsigned short int no_of_nodes_per_element;/**< Number of nodes per element*/
        unsigned short int total_dofs;/**< Total number of degrees of freedom.Derived from number of nodes and the dimension of the problem. Degrees of Freedom are necesarry since we are solving a vector problem and thus at each node we will have a solution in 2 directions or 2 "degrees of freedom"*/
        unsigned short int dofs_per_ele;/**< Number of degrees of freedom per element - Derived from Number of nodes per element and the dimension of the problem*/
        uint8_t quad_rule;/**< Defines the number of Guassian Quadrature points along each direction I use to do the numerical integration. I choose 3 for this problem as this is sufficient to perform the required integration exactly.*/
        uint8_t dim; /**< Defines the dimension of the problem. In this case, it is 2.*/


        std::vector<std::vector<double> > NC; /**<This object is the Nodal Coordinate Matrix. NC[i] gives the x and y - coordinate of the ith global degree of freedom. Note: Each node has 'dimension' degrees of freedom. The size of this matrix will thus be  (No.of nodes, dim)*/
        Eigen::MatrixXd Klocal; /**< Klocal is the elemental stiffness matrix that is used to assemble the Gobal Stiffness. This is a matrix of size (degrees of freedom per element X degrees of freedom per element*/
        std::vector<std::vector<unsigned short int> > EC;/**< Element Connectivity Matrix - EC[i][j] gives the global node number for local node 'j' in element 'i' - Size - (No. of elements, No. of degrees of freedom per element)*/
        std::vector<std::vector<unsigned short int> > EC_2;/**< This is currently only used to write the VTU file. Will be deprecated during submission*/
        std::vector<double> boundary_values; /**<Vector having the dirichlet boundary value wherever its defined and 0 for all other rows - Size - (No. of total degrees of freedom)*/
        std::vector< unsigned short int > boundary_nodes; /**<Vector having all the degree of freedom numbers where there is a dirichlet boundary condition applied - Size depends on the number of dirichlet dofs*/
        std::vector<std::vector<double> > quad_points;/**< Matrix containing the quadrature point locations in terms of \f$(\xi)\f$ and \f$(\eta)\f$. Although there are totally 9 quadrature points, we only define three and loop through these 3 in both directions \f$(\xi)\f$ and \f$(\eta)\f$. Thus, the size of this matrix is (Number of quadrature points in each direction  X Dimension). For more details read source code for function cal_k_local()*/
        std::vector<double> quad_weights;/**<Vector for the weights to be applied at each quadrature point. This is a vector with length = Number of quadrature points in each direction.Again, there should traditionally be 9 weights(total quadrature points), however, we define 3  and loop through them.For more details read source code for function cal_k_local()*/
        Eigen::MatrixXd invJ;/**< In FEM, the basis functions are defined in the parametric space, however, we integrate and find the solution in the real domain. To move from the parametric coordinate system to the real coordinate system, we need the inverse Jacobian.*/
        Eigen::VectorXd U;/**< The solution vector U (displacement).*/
        Eigen::VectorXd F;/**< This is the Forcing Vector that will make up the RHS of the linear equation we finally solve. Since the Dirichlet Conditions are all u = 0 and there is no body force and only a force at a single node, all but one rows of this matrix will be zero.This Vector however is defined as non sparse. This is because, an insertion needs to take place for that one dof where we have a force. Each insertion requires a binary search and it is thus not feaslible to define the F vector a Sparse Vector from the  very beginning. However, before the linear system is solved, this vector will be converted to a Sparse Vector!*/
//        Eigen::SparseMatrix<double> K;
        Eigen::MatrixXd K;/**< This is the Global Stiffness Matrix that makes up the LHS of the linear system of equations we finally solve. By construction, this matrix is largely made up of zeros. Initially it is defined as a non sparse matrix as it requires a lot of data entries and each data entry requires a binary search.However, when the equation is solved for better performance, it is defined as a sparse Matrix. */

};

#endif
