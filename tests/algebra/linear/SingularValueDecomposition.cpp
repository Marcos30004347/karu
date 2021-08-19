// #include <assert.h>

#include "algebra/linear/svd.hpp"
#include "algebra/linear/SingularValueDecomposition.hpp"
// #include <iostream>
int main()
{
	Matrix M(3,3, {
		1.115e-13,       4.677e-06,       -2.342e-03,
		1.119e-05,       -9.764e-05,      3.015e-02,
		-5.584e-03,      7.979e-02,       -3.060e+01,
	});

	Matrix w,s;
	Matrix U,D,V;

	squareSvd(M, 3, 3, w, s);

	// svd(M, U, D, V);
	// printMatrix(U);
	// printMatrix(D);
	// printMatrix(V);
	// dsvd(M, M.rows(), M.columns(), D, V);
	// printMatrix(M);
	// printMatrix(D);
	// printMatrix(V);
// 	double a_[7][9] = {
//     {3.788e+04, 1.288e+05, 1.667e+02, 1.894e+05, 6.439e+05, 8.333e+02, 2.273e+02, 8.333e+02, 1.000e+00},
//     {1.389e+05, 6.944e+05, 8.333e+02, 1.389e+05, 6.944e+05, 8.333e+02, 1.667e+02, 8.333e+02, 1.000e+00},
//     {6.439e+05, 6.439e+05, 7.727e+02, 6.439e+05, 6.439e+05, 7.727e+02, 8.333e+02, 7.727e+02, 1.000e+00},
//     {1.756e+05, 1.756e+05, 2.273e+02, 5.971e+05, 5.971e+05, 7.727e+02, 7.727e+02, 7.727e+02, 1.000e+00},
//     {2.500e+05, 4.000e+05, 5.000e+02, 4.000e+05, 6.400e+05, 8.000e+02, 5.000e+02, 8.000e+02, 1.000e+00},
//     {2.500e+05, 2.939e+05, 5.000e+02, 2.939e+05, 3.456e+05, 5.879e+02, 5.000e+02, 5.879e+02, 1.000e+00},
//     {2.500e+05, 1.879e+05, 5.000e+02, 1.879e+05, 1.412e+05, 3.757e+02, 5.000e+02, 3.757e+02, 1.000e+00},
// 	};

// 	double** a = (double**)malloc(sizeof(double*)*9);
// 		for(int j=0; j<9; j++)
// 			a[j] = (double*)malloc(sizeof(double)*7);
	
// 	for(int i=0; i<7; i++)
// 		for(int j=0; j<9; j++)
// 			a[j][i] = a_[i][j];

// 	double** v = (double**)malloc(sizeof(double*)*9);
// 		for(int j=0; j<9; j++)
// 			v[j] = (double*)malloc(sizeof(double)*9);

// 	double w[9];
	
// 	// svd(a, 7, 9, w, v);

// 	// for(int j=0; j<9; j++)
// 	// {
// 	// 	for(int i=0; i<7; i++)
// 	// 		std::cout << a[j][i] << " ";
// 	// 	std::cout << "\n";
// 	// }
// 	// std::cout << "\n";
// 	// for(int j=0; j<9; j++)
// 	// 		std::cout << w[j] << " ";
// 	// 	std::cout << "\n";



	return 0;
}
