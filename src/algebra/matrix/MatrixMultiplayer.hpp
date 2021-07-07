#pragma once

#include "algebra/matrix/MatrixData.hpp"

namespace karu {
namespace algebra {

class MatrixMultiplayer {
	public:
	static MatrixData multiply(MatrixData* A, MatrixData* B, bool A_T, bool B_T);
};

}
}
