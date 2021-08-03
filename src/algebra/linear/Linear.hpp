#pragma once

#include "algebra/core/types.hpp"
#include "algebra/matrix/Matrix.hpp"
#include <cmath>

namespace karu::algebra {

#define pi 3.141592653589793238462643383279502884197169399375105820974944

Matrix cross(Matrix A, Matrix B);
f32 norm(Matrix v);
f32 radians(f32 degrees);
f32 roundToPrecision(f32 val, u64 p);

}
