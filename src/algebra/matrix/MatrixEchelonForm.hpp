#pragma once

#include "algebra/matrix/MatrixData.hpp"

namespace karu::algebra {

class MatrixEchelonForm {

public:
	static void toEchelonForm(MatrixData* a);
};

}
