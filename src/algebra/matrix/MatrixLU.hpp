#pragma once

#include "algebra/matrix/MatrixData.hpp"

namespace karu::algebra {

class MatrixLU {
	public:
	static void LUdecompose(MatrixData* L, MatrixData* U, const MatrixData* const A);
	static i32 LUPdecompose(MatrixData* A, MatrixData* P);
	static i32 LUPSolve(const MatrixData* const A, const MatrixData* const P,  const MatrixData* const b,  MatrixData* x);
	static i32 LUPInvet(const MatrixData* const A, const MatrixData* const P,  MatrixData* A_Inv);
	static f32 LUPDeterminant(const MatrixData* const A, const MatrixData* const P);
};

}
