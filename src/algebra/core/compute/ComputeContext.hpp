#pragma once

#include "algebra/core/types.hpp"
#include "algebra/core/compute/OpenCL.hpp"

namespace karu {
namespace algebra {
namespace compute {

class Context {
	cl_command_queue queue;
	cl_context context;

	public:

	Context();
	~Context();

	cl_context getComputeContext();
	cl_command_queue getComputeQueue();

	static Context* initContext();
};

}
}
}
