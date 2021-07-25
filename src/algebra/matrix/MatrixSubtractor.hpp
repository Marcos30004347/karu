#pragma once

#include "algebra/matrix/MatrixData.hpp"

namespace karu {
namespace algebra {

class MatrixSubtractor {
	public:
	static void sub(MatrixData* C, const MatrixData* const A, const MatrixData* const B, bool A_T, bool B_T);
};

}
}
