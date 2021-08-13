#include <assert.h>
#include <iostream>
#include <chrono>

#include "algebra/matrix/Matrix.hpp"
#include "algebra/linear/SingularValueDecomposition.hpp"

using namespace karu::algebra;


void echelonFormTests()
{
	Matrix A = Matrix(3, 4, {
		1, 2, -1, -4,
		2, 3, -1, -11,
		-2, 0, -3, 22
	});
	Matrix A_echelon = echelonForm(A);

	assert(A_echelon[0][0] == 1.f);
	assert(A_echelon[1][0] == 0.f);
	assert(A_echelon[2][0] == 0.f);

	assert(A_echelon[0][1] == 0.f);
	assert(A_echelon[1][1] == 1.f);
	assert(A_echelon[2][1] == 0.f);

	assert(A_echelon[0][2] == 0.f);
	assert(A_echelon[1][2] == 0.f);
	assert(A_echelon[2][2] == 1.f);

	assert(A_echelon[0][3] == -8.f);
	assert(A_echelon[1][3] == 1.f);
	assert(A_echelon[2][3] == -2.f);
}

void nullspaceTests()
{
	Matrix A = Matrix(3, 4, {
		1, 2, 3, 4,
		1, 3, 5, 6,
		2, 5, 8, 10
	});

	Matrix A_space = nullspace(A);

	assert(A_space[0][0] == 1);
	assert(A_space[0][1] == -2);
	assert(A_space[0][2] == 1);
	assert(A_space[0][3] == 0);

	assert(A_space[1][0] == 0);
	assert(A_space[1][1] == -2);
	assert(A_space[1][2] == 0);
	assert(A_space[1][3] == 1);

	Matrix B = Matrix(2, 4, {
		-1, 1, 2, 4,
		2, 0, 1, -7,
	});

	Matrix B_space = nullspace(B);

	assert(B_space[0][0] == -0.5);
	assert(B_space[0][1] == -2.5);
	assert(B_space[0][2] == 1);
	assert(B_space[0][3] == 0);

	assert(B_space[1][0] == 3.5);
	assert(B_space[1][1] == -0.5);
	assert(B_space[1][2] == 0);
	assert(B_space[1][3] == 1);

	Matrix C = Matrix(2, 2, {
		2, 1,
		1, 2,
	});

	Matrix C_space = nullspace(C);

	assert(C_space[0][0] == 0);
	assert(C_space[0][1] == 0);
}

// void svdTests()
// {
// 	std::chrono::steady_clock::time_point begin;
// 	std::chrono::steady_clock::time_point end;

// 	float* w = vector(1, 3);

// 	float** m = matrix(1, 2, 1, 3);
// 	float** v = matrix(1, 3, 1, 3);
	
// 	m[1][1] = 3; m[1][2] = 2; m[1][3] = 2;
// 	m[2][1] = 2; m[2][2] = 3; m[2][3] = -2;
	
// 	std::cout << m[1][1] << " "; std::cout << m[1][2] << " "; std::cout << m[1][3] << "\n";
// 	std::cout << m[2][1] << " "; std::cout << m[2][2] << " "; std::cout << m[2][3] << "\n";

// 	begin = std::chrono::steady_clock::now();
// 	svdcmp(m, 2, 3, w, v);
// 	end = std::chrono::steady_clock::now();

// 	std::cout << "Naive Total: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "µs" << std::endl;
// 	std::cout << "Naive Total: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms" << std::endl;

// 	std::cout << "\n";
// 	std::cout << m[1][1] << " "; std::cout << m[1][2] << " "; std::cout << m[1][3] << "\n";
// 	std::cout << m[2][1] << " "; std::cout << m[2][2] << " "; std::cout << m[2][3] << "\n";
// 	std::cout << "\n";
	
// 	std::cout << w[1] << " "; std::cout << 0 << " "; std::cout <<  0 << "\n";
// 	std::cout << 0 << " "; std::cout << w[2] << " "; std::cout <<  0 << "\n";
// 	std::cout << 0 << " "; std::cout << 0 << " "; std::cout << w[3] << "\n";
// 	std::cout << "\n";

// 	std::cout << v[1][1] << " "; std::cout << v[2][1] << " "; std::cout << v[3][1] << "\n";
// 	std::cout << v[1][2] << " "; std::cout << v[2][2] << " "; std::cout << v[3][2] << "\n";
// 	std::cout << v[1][3] << " "; std::cout << v[2][3] << " "; std::cout << v[3][3] << "\n";
// 	std::cout << "\n";


// 	Matrix U(2,3, {
// 		m[1][1], m[1][2], m[1][3],
// 		m[2][1], m[2][2], m[2][3],
// 	});

// 	Matrix W(3,3, {
// 		w[1], 0, 0,
// 		0, w[2], 0,
// 		0, 0, w[3],
// 	});

// 	Matrix V(3,3, {
// 		v[1][1], v[1][2], v[1][3],
// 		v[2][1], v[2][2], v[2][3],
// 		v[3][1], v[3][2], v[3][3],
// 	});

// 	Matrix A = U*W*transpose(V);
// 	printMatrix(A);
// }



void _svdTests()
{
	Matrix M(2,3, {
		3, 2, 2,
		2, 3, -2
	});
	
	Matrix U(2,3);
	Matrix V(3,3);
	Matrix D(3,1);
	
	svd(M, U, D, V);

	Matrix W(3,3, {
		D[0][0], 0, 0,
		0, D[1][0], 0,
		0, 0, D[2][0],
	});

	Matrix A = U*W*transpose(V);
	

	assert(fabs(A[0][0] - M[0][0]) < 0.000009);
	assert(fabs(A[0][1] - M[0][1]) < 0.000009);
	assert(fabs(A[0][2] - M[0][2]) < 0.000009);
	assert(fabs(A[1][0] - M[1][0]) < 0.000009);
	assert(fabs(A[1][1] - M[1][1]) < 0.000009);
	assert(fabs(A[1][2] - M[1][2]) < 0.000009);


}


int main()
{
	Matrix A(3, 3, {1, 2, 3, 4, 5, 6, 7, 8, 9}, 2, 2);
	
	int v = 1;
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
			assert(A[i][j] == v++);

	Matrix B(9, 9, {
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
	}, 4, 4);

	for(int i=0; i<9; i++)
		for(int j=0; j<9; j++)
			assert(B[i][j] == j+1);


	Matrix C(9, 9, {
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
	}, 4, 4);

	Matrix D(9, 9, {
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
	}, 4, 4);

	Matrix E = C + D;

	for(int i=0; i<9; i++)
		for(int j=0; j<9; j++)
			assert(E[i][j] == C[i][j] + D[i][j]);

	Matrix F(9, 9, {
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
	}, 2, 2);

	Matrix G = C + F;

	for(int i=0; i<9; i++)
		for(int j=0; j<9; j++)
			assert(G[i][j] == C[i][j] + F[i][j]);

	Matrix H = C*D;

	float data[] = {45, 90, 135, 180, 225, 270, 315, 360, 405};

	for(int i=0; i<9; i++)
		for(int j=0; j<9; j++)
			assert(H[i][j] == data[j]);

	Matrix I = C*F;

	for(int i=0; i<9; i++)
		for(int j=0; j<9; j++)
			assert(I[i][j] == data[j]);


	Matrix J(9, 9, {
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
		1, 2, 3, 4, 5, 6, 7, 8, 9,
	}, 4, 4);

	Matrix K(9, 1, {
		1, 2, 3, 4, 5, 6, 7, 8 ,9
	}, 4, 1);

	Matrix L = J*K;

	assert(L.rows() == 9);
	assert(L.columns() == 1);

	for(int i=0; i<9; i++)
		assert(L[i][0] == 285);

	Matrix M = transpose(J);

	assert(M.rows() == J.columns());
	assert(M.columns() == J.rows());

	for(int i=0; i<M.rows(); i++)
		for(int j=0; j<M.columns(); j++)
			assert(M[i][j] == J[j][i]);

	Matrix N = transpose(L);
	assert(N.rows() == L.columns());
	assert(N.columns() == L.rows());

	for(int i=0; i<N.rows(); i++)
		for(int j=0; j<N.columns(); j++)
			assert(N[i][j] == L[j][i]);

	Matrix P = N*L;

	assert(P.columns() == 1);
	assert(P.rows() == 1);
	assert(P[0][0] == 731025);

	Matrix Q = M*3.f;

	for(int i=0; i<Q.rows(); i++)
		for(int j=0; j<Q.columns(); j++)
			assert(Q[i][j] == 3*M[i][j]);

	Matrix W = Q/3.f;

	for(int i=0; i<Q.rows(); i++)
		for(int j=0; j<Q.columns(); j++)
			assert(W[i][j] == Q[i][j]/3.f);
	
	Matrix X = 3;

	Matrix Z = Q/X;

	for(int i=0; i<Q.rows(); i++)
		for(int j=0; j<Q.columns(); j++)
			assert(Z[i][j] == Q[i][j]/3.f);
	
	Matrix Y = W*X;

	for(int i=0; i<Q.rows(); i++)
		for(int j=0; j<Q.columns(); j++)
			assert(Y[i][j] == W[i][j]*3.f);
	
	Matrix T(3,3,{
		1,1,0,
		2,1,3,
		3,1,1
	});

	std::pair<Matrix, Matrix> T_LU = LUDecomposition(&T);

	assert(T_LU.first[0][0] == 1);
	assert(T_LU.first[0][1] == 0);
	assert(T_LU.first[0][2] == 0);
	assert(T_LU.first[1][0] == 2);
	assert(T_LU.first[1][1] == -1);
	assert(T_LU.first[1][2] == 0);
	assert(T_LU.first[2][0] == 3);
	assert(T_LU.first[2][1] == -2);
	assert(T_LU.first[2][2] == -5);

	assert(T_LU.second[0][0] == 1);
	assert(T_LU.second[0][1] == 1);
	assert(T_LU.second[0][2] == 0);
	assert(T_LU.second[1][0] == 0);
	assert(T_LU.second[1][1] == 1);
	assert(T_LU.second[1][2] == -3);
	assert(T_LU.second[2][0] == 0);
	assert(T_LU.second[2][1] == 0);
	assert(T_LU.second[2][2] == 1);

	Matrix System(2,2,{
		3,2,
		2,6,
	});

	// LU decomposition tests
	std::pair<Matrix, Matrix> System_Permutation = LUPDecomposition(System);
	Matrix b(2, 1, {2,-8});

	Matrix x = LUPSolve(System_Permutation.first, System_Permutation.second, b);

	assert(x[0][0] == 2);
	assert(x[1][0] == -2);

	Matrix System_Inv = LUPInverse(System_Permutation.first, System_Permutation.second);
	Matrix Res = System*System_Inv;
	assert(Res[0][0] - 1 < 0.00000001);
	assert(Res[0][1] - 0 < 0.00000001);
	assert(Res[1][0] - 0 < 0.00000001);
	assert(Res[1][1] - 1 < 0.00000001);

	karu::f32 System_det = LUPDeterminant(System_Permutation.first, System_Permutation.second);
	assert(System_det == 14);
	

	echelonFormTests();
	nullspaceTests();
	_svdTests();
	return 0;
}
