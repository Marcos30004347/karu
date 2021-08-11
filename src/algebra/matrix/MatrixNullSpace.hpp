#pragma once

#include "algebra/matrix/MatrixData.hpp"

namespace karu::algebra {

class MatrixNullSpace 
{
public:
	static void nullSpace(MatrixData* M, MatrixData& ns);
};

}
