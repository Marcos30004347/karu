#pragma once

#include "algebra/core/types.hpp"
#include "algebra/core/compute/OpenCL.hpp"
#include "algebra/core/compute/ComputeBuffer.hpp"

namespace karu {
namespace algebra {
namespace compute {

class Program {
	public:

	Program(const char* src);
	~Program();

	cl_program getProgram();
	
	private:

	const char* cp_src;
	cl_program cp_program;	
};

}
}
}

