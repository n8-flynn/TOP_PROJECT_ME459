#include "top.h"
#include <sstream>

using namespace Eigen;

using namespace std; 

/*! \mainpage Main functions
*	 main.cpp
* \n top.cpp
* \n FE.cpp
* \n OC.cpp
* \n check.cpp
*/

void writeToCsv(string fileName, MatrixXd  matrix)
{
	//! https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html

	const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");
	ofstream file(fileName);
	if (file.is_open())
	{
		file << matrix.format(CSVFormat);
		file.close();
	}
}

MatrixXd mfilter(MatrixXd &m1, double filter) {
	int r = m1.rows();
	int c = m1.cols();


	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			if (m1(i, j) > filter) {
				m1(i, j) = 1;
			}
			if (m1(i, j) < filter) {
				m1(i, j) = 0;
			}
		}
	}
	return m1;
}

int main(int argc, char* argv[]) 
{
    if(argc < 6){
        cout<<"Please provide sufficient inputs to run the topology code"<<endl;
        return 0;
    }
    istringstream x(argv[1]);
    size_t nelx;
    x >> nelx;
    
    istringstream y(argv[2]);
    size_t nely;
    y >> nely;
    double volfrac = stod(argv[3]);
    double penal = stod(argv[4]);
    double rmin = stod(argv[5]);
    
//    istringstream wh_(argv[6]);
    int wh = atoi(argv[6]);
//    wh_ >> wh;

    if(wh > 2){
        cout<<"Please enter boundary condition number from 0 to 2. Anything else is invalid";
        return 0;
    }
    cout<<"You have provided nelx = "<<nelx<<", nely = "<<nely<<", volfrac = "<<volfrac<<", penal = "<<penal<<", rmin = "<<rmin<<" and you have choosen the boundary condition of type "<<wh<<endl;
	
	MatrixXd output = top(nelx, nely, volfrac, penal, rmin,wh);

	//mfilter(output, 0.4);

	//cout << output << endl;
	
	writeToCsv("density_field.csv", output);
	
	return 0;
}
