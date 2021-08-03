#include "algebra/linear/Linear.hpp"




namespace karu::algebra {

Matrix cross(Matrix A, Matrix B)
{
	return Matrix(3,1,{
		A[1][0]*B[2][0] - A[2][0]*B[1][0],
		A[2][0]*B[0][0] - A[0][0]*B[2][0],
		A[0][0]*B[1][0] - A[1][0]*B[0][0],
	});
}

f32 norm(Matrix v) {
	return std::sqrt( v[0][0]*v[0][0] + v[1][0]*v[1][0] + v[2][0]*v[2][0] ) ;
}

f32 radians(f32 degrees)
{
	return ( degrees * pi) / 180.0 ;
}


f32 roundToPrecision(f32 val, u64 p)
{
	f32 k = 1;
	while(p--) k*=10;

	val *= k;
	return ((long long)val)/k;
}

}
