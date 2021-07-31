#pragma once

#include "algebra/core/types.hpp"
#include "algebra/matrix/Matrix.hpp"
#include <cmath>

namespace karu::algebra {

Matrix cross(Matrix A, Matrix B)
{
	return {
		A[1][0]*B[2][0] - A[2][0]*B[1][0],
		A[2][0]*B[0][0] - A[0][0]*B[2][0],
		A[0][0]*B[1][0] - A[1][0]*B[0][0],
	};
}

f32 norm(Matrix v) {
	return std::sqrt( v[0][0]*v[0][0] + v[1][0]*v[1][0] + v[2][0]*v[2][0] ) ;
}


}