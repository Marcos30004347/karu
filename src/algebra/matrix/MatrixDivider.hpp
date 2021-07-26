#pragma once

#include "algebra/matrix/MatrixData.hpp"

namespace karu {
namespace algebra {

class MatrixDivider {
	public:
	static void div(MatrixData* C, const MatrixData* const A, const MatrixData* const B, bool A_T, bool B_T);
	static void div(MatrixData* C, const MatrixData* const A, const float alpha);
};

}
}
