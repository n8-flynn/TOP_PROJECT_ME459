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
	size_t nelx = 32;
	size_t nely = 20;
	double volfrac = 0.4;
	double penal = 3;
	double rmin = 1.2;
	
	printf("Top starting\n");
	
	MatrixXd output = top(nelx, nely, volfrac, penal, rmin);

	//mfilter(output, 0.5);

	cout << output << endl;
	
	writeToCsv("density_field.csv", output);
	
	return 0;
}
