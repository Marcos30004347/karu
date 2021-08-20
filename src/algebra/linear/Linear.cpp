#include "algebra/linear/Linear.hpp"

namespace karu::algebra {

Matrix outerProduct(Matrix& v0, Matrix& v1)
{
	Matrix o(v0.rows(), v1.rows());

	for(i32 i=0; i<v0.rows(); i++)
	{
		for(i32 j=0; j<v1.rows(); j++)
		{
			o[i][j] = v0[i][0] * v1[j][0];
		}
	}
	return o;
}	

Matrix cross(Matrix A, Matrix B)
{
	return Matrix(3,1,{
		A[1][0]*B[2][0] - A[2][0]*B[1][0],
		A[2][0]*B[0][0] - A[0][0]*B[2][0],
		A[0][0]*B[1][0] - A[1][0]*B[0][0],
	});
}

f32 norm(Matrix& A){	
	f32 n = 0;
	
	for(i32 l=0; l<A.rows(); l++)
		n += A[l][0]*A[l][0];
	
	return std::sqrt(n);
}

f32 norm(Matrix A){	
	f32 n = 0;
	
	for(i32 l=0; l<A.rows(); l++)
		n += A[l][0]*A[l][0];
	
	return std::sqrt(n);
}
// f32 norm(Matrix v) {
// 	return std::sqrt( v[0][0]*v[0][0] + v[1][0]*v[1][0] + v[2][0]*v[2][0] ) ;
// }

f32 radians(f32 degrees)
{
	return ( degrees * PI) / 180.0 ;
}


f32 roundToPrecision(f32 val, u64 p)
{
	f32 k = 1;
	while(p--) k*=10;

	val *= k;
	return ((long long)val)/k;
}

}
