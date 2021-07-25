#pragma once

#include "algebra/matrix/MatrixData.hpp"

namespace karu {
namespace algebra {

class MatrixAdder {
	public:
	static void add(MatrixData* C, const MatrixData* const A, const MatrixData* const B, bool A_T, bool B_T);
};

// convulação
// recude


}
}
