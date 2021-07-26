#pragma once

#include "algebra/matrix/MatrixData.hpp"

namespace karu::algebra {

class MatrixTransposer {
	public:
	static void transpose(MatrixData* C, const MatrixData* const A);
};

}
