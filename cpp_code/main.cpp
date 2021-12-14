#include "top.h"

using namespace Eigen;

using namespace std; 

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

int main(int argc, char* argv[]) 
{
	size_t nelx = 20;
	size_t nely = 10;
	double volfrac = 0.5;
	double penal = 1.5;
	double rmin = 1.125;
	
	printf("Top starting\n");
	
	MatrixXd output = top(nelx, nely, volfrac, penal, rmin);
	
	printf("Top done\n");
	
	writeToCsv("density_field.csv", output);
	return 0;
}
