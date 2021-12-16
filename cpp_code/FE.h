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
        FE(unsigned short int nelx,unsigned short int nely,double length, double breadth,double youngs_mod, double pois_rat); /**< The constructor takes the number of elements along the x axis (nelx), the number of elements along the Y axis (nely), the length and breadth (although these are set to always equal nelx and nely), the Youngs Modulus of the material (youngs_mod) and the Poisson Ratio (pois_rat).All of these parameters are then assigned to their respective class members.*/
        double basis_function(unsigned short int node, double xi,double eta); /**<Calculates the value basis function corresponding to the node "node" and at the point 'xi' and 'eta' in the parametric space.'xi' and 'eta' are the coordinates of the parametric space. This method is used while finding the Elemental stiffness matrix (klocal).*/
        std::vector<double> basis_gradient(unsigned short int node, double xi,double eta); /**<Calculates the gradient of the basis function corresponding to the node "node" and at the point 'xi' and 'eta' with respect to both 'xi' and 'eta'. Thus, this method returns a vector of length 2 whose 0th element is the gradient with respect to 'xi' and the 1st element is the gradient with respect to 'eta'.*/
        void mesh(uint8_t no_quad_points); /**< This method fills up the Nodal Coordinate matrix and the Elemental Connectivity matrix whihc define the mesh of the domain. It also takes as input the number of quadrature points in order to define the location of the quadrature points and the quadrature weights. Currently this method can only take 3 quadrature points, however, this can easily be explanded in the future.*/
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
        unsigned short int nelx_,nely_,nel,nnx_,nny_,no_of_nodes,no_of_nodes_per_element,total_dofs,dofs_per_ele;
        uint8_t quad_rule,dim; // standard mesh descriptions


        std::vector<std::vector<double> > NC; //Nodal Connectivity - NC[i] gives the x and y - coordinate of the ith global node. Size - (No.of nodes, dim)
//        std::vector<std::vector<double> > Klocal;
        Eigen::MatrixXd Klocal;
        std::vector<std::vector<unsigned short int> > EC;
        std::vector<std::vector<unsigned short int> > EC_2;//Elemental connectivity - EC[i][j] gives the global node number for local node j in element i - Size - (No. of elements, No. of nodes per element)

        std::vector<double> boundary_values; // Vector having the dirichlet boundary value wherever its defined and 0 in all other entries - Size (No. of nodes)
        std::vector< unsigned short int > boundary_nodes; //Vector having all the nodes that are part of the dirichlet boundary - Size depends on the number of dirichlet nodes
        std::vector<std::vector<double> > quad_points; // Vector for the location of the quad points in terms of xi
        std::vector<double> quad_weights; // Vector for the weights at the quad points. Still 1D as the weigths do not depend on xi or eta
        Eigen::MatrixXd invJ; // Inverse Jacobian needed
        Eigen::VectorXd U; // Using Eigne vector to define to solution for ease of solving the linear equation
        Eigen::VectorXd F; // No forcing but the dirichlet conditions will apply
//        Eigen::MatrixXd K; // Not using sparse matrix to get maximum speed irrespective of memory usage
        Eigen::SparseMatrix<double> K; // The global stiffness matrix - Sparse because its a tridiagonal matrix and so alot of elements are 0 - saves memory

};

#endif
