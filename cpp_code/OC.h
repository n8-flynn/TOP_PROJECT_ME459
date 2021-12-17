#ifndef OC_h
#define OC_h

#include <iostream> 
#include <algorithm>
#include <cmath>
#include <stddef.h>
#include "Eigen"

using namespace Eigen;
using namespace std;

MatrixXd OC(unsigned int nelx, unsigned int nely, double volfrac, MatrixXd x, MatrixXd dc, double move_val);

/*!
	\fn MatrixXd OC(unsigned int nelx, unsigned int nely, double volfrac, MatrixXd x, MatrixXd dc, double move_val)
	\brief Optimization function. 
	\param nelx Number of elements in the x direction.
	\param nely Number of elements in the y direction.	
	\param volfrac Volume fraction.
	\param x Volume fraction field that is nely by nelx. 
	\param dc (not sure on this).
	\param move_val Distance volume fraction matrix is indexed (not sure on this).
*/

MatrixXd mmax(MatrixXd m1, MatrixXd m2, MatrixXd m3);

/*!
	\fn MatrixXd mmax(MatrixXd m1, MatrixXd m2, MatrixXd m3);
	\brief Function used to find the max value between two input arrays, element by element. Returns array of max values between both arrays. 
	\param m1 Input array 1. 
	\param m2 Input array 2.
	\param m3 Output array. 
*/

MatrixXd mmin(MatrixXd m1, MatrixXd m2, MatrixXd m3);

/*!
	\fn MatrixXd mmax(MatrixXd m1, MatrixXd m2, MatrixXd m3);
	\brief Function used to find the max value between two input arrays, element by element. Returns array of max values between both arrays. 
	\param m1 Input array 1. 
	\param m2 Input array 2.
	\param m3 Output array.
*/

void writeToCsv2(string fileName, MatrixXd  matrix);

/*!
	\fn void writeToCsv2(string fileName, MatrixXd  matrix);
	\brief Writes inputed eigen arrrays to csv files. 
	\param filename Desired file name. 
	\param matrix Inputed eigen array.
*/

#endif
