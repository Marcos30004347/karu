#include <cmath>
#include <assert.h>
#include "algebra/matrix/Matrix.hpp"

namespace karu::algebra 
{

float norm(Matrix& A){
	assert(A.columns() == 1);
	
	f32 n = 0;
	
	for(i32 l=0; l<A.rows(); l++)
		n += A[l][0]*A[l][0];
	
	return std::sqrt(n);
}

Matrix _conjGrad(Matrix A, Matrix x, Matrix b, float tolerance)
{
	Matrix r = b - A*x;
	Matrix p = r;

	Matrix rsold = transpose(r)*r;

	while(norm(r) > tolerance)
	{
		std::cout << norm(r) << std::endl;
	
		Matrix Ap = A*p;
		Matrix alpha = rsold/(transpose(p)*Ap);
	
		x = x + alpha*p;
		r = r - alpha*Ap;
		
		Matrix rsnew = transpose(r)*r;
	
		std::cout << norm(r) << std::endl;
	
		// if(norm(r) < tolerance)
		// {
		// 	break;
		// }

		p = r + (rsnew/rsold)*p;
		rsold = rsnew;
	}
	return x;
}


// Solve Symmetric Positive Definite systems A*x = b, returns x
Matrix conjGrad(Matrix A, Matrix x, Matrix b, float tolerance)
{
	Matrix p;
	Matrix r[2];

	r[0] = b - A*x;
	p = r[0];

	Matrix W = A*p;
 	Matrix alpha = (transpose(r[0])*r[0])/(transpose(p)*W);

	x = x + alpha*p;
	r[1] = r[0] - alpha*W;

	i32 k = 1;

	while (norm(r[1]) > tolerance)
	{
		Matrix beta = (transpose(r[1])*r[1])/(transpose(r[0])*r[0]);
	
		p = r[1] + beta*p;

		W = A*p;

		alpha =(transpose(r[1])*r[1])/(transpose(p)*W);

		x = x + alpha*p;

		r[0] = r[1];
		r[1] = r[1] - alpha*W;

		k = k+1;	
	}

	return x;
}



// Solve Symmetric Positive Definite systems A*x = b, returns x
Matrix preConjGrad(Matrix A, Matrix x, Matrix b, Matrix M, float tolerance)
{
	Matrix p;

	Matrix r[2];
	Matrix z[2];

	r[0] = b - A*x;
 	
	std::pair<Matrix, Matrix> LUP = LUPDecomposition(M);
	
	z[0] = LUPSolve(LUP.first, LUP.second, r[0]);

	p = z[0];

	Matrix W = A*p;

 	Matrix alpha = (transpose(r[0])*z[0])/(transpose(p)*W);

	x = x + alpha*p;
	r[1] = r[0] - alpha*W;

	i32 k = 1;

	while (norm(r[1]) > tolerance)
	{
		z[1] = LUPSolve(LUP.first, LUP.second, r[1]);

		Matrix beta = (transpose(r[1])*z[1])/(transpose(r[0])*z[0]);
	
		p = z[1] + beta*p;

		W = A*p;

		alpha =(transpose(r[1])*z[1])/(transpose(p)*W);

		x = x + alpha*p;

		z[0] = z[1];
		r[0] = r[1];

		r[1] = r[1] - alpha*W;

		k = k+1;	
	}

	return x;
}

}
