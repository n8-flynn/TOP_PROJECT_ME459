// Created by Huzaifa Mustafa Unjhawala

#ifndef FE_H
#define FE_H

// All the required libraries
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>



typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
class FE{
	FEM(unsigned int nelx,unsigned int nely,unsigned int length, unsigned int breadth,double penal,double youngs_mod, double pois_rat); // The constructor takes the required arguments
	double basis_function(unsigned int node, double xi,double eta); //Calculates the basis function corresponding to the node "node" and at the point 'xi' and 'eta' in the parametric space
	std::vector<double> basis_gradient(unsigned int node, double xi,double eta); //Calculates the gradient of the basis function - similar to above 
	void mesh(unsigned int elements_l,unsigned int elements_b,double length, double breadth,unsigned int no_quad_points); // Function to mesh the domain - Fills the Nodal connectivity and the Elemental Conncectivity matrices - Can even handle different number of elements along each axis
	// As of now, user will only have ability to define where the boundary conditions act from outside the code but can change the values of boundary conditions from outside the code
	void define_boundary_condition(double force, double g); // Function to fill the boundary_values (stores the values at the boundaries) and the boundary_nodes (stores the global node number of the nodes on the boundary)
	double C(unsigned int i, unsigned int j, unsigned int k, unsigned int l); // Used to get the elasticity tensor C
	void init_data_structs(); //To resize all the global matrices based on the mesh - Internal function
	void fe_impl(Eigen::MatrixXd x); //Does the local looping , the assembly and applying the boundary conditions - fills K,F and M which is then used in the steady and transient state solution. Also applies the dirichlet conditions
	Eigen::VectorXd solve(); // Solves and returns U which is then used in the toplogy code


	// Class datastructures
	double L,B,g1,E,nu,lambda,mu; //Standard constants - L - Length, B - breadth, g1 - Dirichlet conditon, E - Youngs Modulus, nu - Poissons ration, lambda and mu
	// are the lame's parameters 
	unsigned int nelx_,nely_,nel,nnx,nny,nn,no_of_nodes_per_element,quad_rule; // standard mesh descriptions


    std::vector<std::vector<double> > NC; //Nodal Connectivity - NC[i] gives the x and y - coordinate of the ith global node. Size - (No.of nodes, dim)
	std::vector<std::vector<int> > EC; //Elemental connectivity - EC[i][j] gives the global node number for local node j in element i - Size - (No. of elements, No. of nodes per element)
	std::vector<double> boundary_values; // Vector having the dirichlet boundary value wherever its defined and 0 in all other entries - Size (No. of nodes)
    std::vector< unsigned int > boundary_nodes; //Vector having all the nodes that are part of the dirichlet boundary - Size depends on the number of dirichlet nodes
	std::vector<double> quad_points; // Vector for the location of the quad points in terms of xi
	std::vector<double> quad_weights; // Vector for the weights at the quad points. Still 1D as the weigths do not depend on xi or eta
    Eigen::VectorXd U; // Using Eigne vector to define to solution for ease of solving the linear equation
    Eigen::VectorXd F; // No forcing but the dirichlet conditions will apply
    Eigen::MatrixXd K; // Not using sparse matrix to get maximum speed irrespective of memory usage

};


