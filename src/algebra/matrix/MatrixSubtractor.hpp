#pragma once

#include "algebra/matrix/MatrixData.hpp"

namespace karu {
namespace algebra {

class MatrixSubtractor {
	public:
	static MatrixData sub(MatrixData* A, MatrixData* B, bool A_T, bool B_T);
};

}
}
