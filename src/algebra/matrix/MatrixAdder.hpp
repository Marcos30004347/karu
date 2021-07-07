#pragma once

#include "algebra/matrix/MatrixData.hpp"

namespace karu {
namespace algebra {

class MatrixAdder {
	public:
	static MatrixData add(MatrixData* A, MatrixData* B, bool A_T, bool B_T);
};

}
}
