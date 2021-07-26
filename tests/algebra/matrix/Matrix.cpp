#include <assert.h>
#include <iostream>

#include "algebra/matrix/Matrix.hpp"

using namespace karu::algebra;

void print(Matrix* A)
{
	for(int i=0;i<A->m_data.m_stored_lines; i++)
	{
		for(int j=0; j<A->m_data.m_stored_column;j++)
		{
			std::cout << A->m_data.get(i,j) << " ";
		}
		std::cout << std::endl;
	}
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
	
	return 0;
}
