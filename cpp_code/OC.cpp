#include "OC.h"
#include <fstream>
#include <iostream>

using namespace Eigen;
using namespace std;

void writeToCsv2(string fileName, MatrixXd  matrix) {
    /*! 
        \brief Writes inputed eigen arrrays to csv files.
        \n Source: https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html 
        \param filename Desired file name.
        \param matrix Inputed eigen array.
    */
    
    const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");
    ofstream file(fileName);
    if (file.is_open())
    {
        file << matrix.format(CSVFormat);
        file.close();
    }
}

MatrixXd mmin(MatrixXd m1, MatrixXd m2,MatrixXd m3) {
    /*!
        \brief Function used to find the min value between two input arrays at the same location, element by element. Returns array of min values between both arrays.
        \param m1 Input array 1.
        \param m2 Input array 2.
        \param m3 Output array.
    */

    size_t r = m1.rows();
    size_t c = m1.cols();

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            if (m1(i, j) < m2(i, j)) {
                m3(i,j) = m1(i,j);
            }
            else {
                m3(i, j) = m2(i, j);
            }
        }
    }
    return m3;
}

MatrixXd mmax(MatrixXd m1, MatrixXd m2, MatrixXd m3) {
    /*!
        \brief Function used to find the max value between two input arrays, element by element. Returns array of max values between both arrays.
        \param m1 Input array 1.
        \param m2 Input array 2.
        \param m3 Output array.
    */
    
    size_t r = m1.rows();
    size_t c = m1.cols();
    
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            if (m1(i, j) > m2(i, j)) {
                m3(i,j) = m1(i,j);
            }
            else {
				m3(i,j) = m2(i,j);
            }
        }
    }
    return m3;
}

MatrixXd OC(unsigned int nelx, unsigned int nely, double volfrac,MatrixXd x,MatrixXd dc, double move_val) {
    /*!
        \brief Optimization criteria based optimizer. Based off of Solid Isotropic Material with Penalization (SIMP) method developed by Bendsoe and Sigmund. 
        *Using this function will predict the best material distribution within the x Matrix. 
        *Assigns densities of 1 and 0 to the domain of x.
        *Material is needed at densities equal to 1 and material is not needed at densities of 0.
        \param nelx Number of elements in the x direction.
        \param nely Number of elements in the y direction.
        \param volfrac Volume fraction.
        \param x Volume fraction field that is nely by nelx.
        \param dc Change in complaince with respect to density. Gives the sensitivity of the compliance.
        \param move_val Move limit (move_val) is nothing but the maximum fractional change allowd for each topology design. 
    */
    
    double l1 = 0.0; //Lower bound of the lagrangian multiplier. 
	double l2 = 100000; //Upper boud of the lagrandgian multiplier. 
	double lmid; //Bound half way beteen the lower and upper bound. 
	double op1 = 0.001; //Optimal criteria #1
	double op2 = 1.0; //Optimal criteria #2
	
	MatrixXd move(nely, nelx); 
    move.setConstant(move_val);
	
    MatrixXd xnew(nely, nelx);
	MatrixXd newdc(nely, nelx);
	
    MatrixXd op1m(nely, nelx); 
    op1m.setConstant(op1);

    MatrixXd op2m(nely, nelx); 
    op2m.setConstant(op2);

	while (l2 - l1 > 0.0001) { //Convergence criteria.
		lmid = 0.5 * (l2 + l1);
       
        MatrixXd r1 = x.array()*((-dc/lmid).cwiseSqrt()).array();
        MatrixXd r2 = x + move;
        
        MatrixXd r3(nely,nelx); 
        r3 = mmin(r2,r1,r3); 
        
        MatrixXd r4(nely,nelx); 
        r4 = mmin(r3,op2m,r4);
        
        MatrixXd r5 = x - move;
        
        MatrixXd r6(nely,nelx);
        r6 = mmax(r5,r4,r6);
        
        MatrixXd r7(nely,nelx);
        r7 = mmax(op1m,r6,r7);
        xnew = r7;
        
		if (xnew.sum() - volfrac * nelx * nely > 0) {
            //Used to increase or degree the bounds of the lagrangian multiplier.
            //If the sum of the densities in x is larger then the total density of the domain, the bounds become between lmid and l2.
            //If the sume of the densities in x is less then the total density of the domain, the bounds become between l1 and lmid. 
			l1 = lmid;
		}
		else
			l2 = lmid;
	}
	return xnew; //Returns the optimized / updated density field x as xnew. 
}


