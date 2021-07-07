#pragma once

#include <vector>

#include "algebra/core/types.hpp"
#include "algebra/core/compute/OpenCL.hpp"

namespace karu {
namespace algebra {
namespace compute {

class Context {
	std::vector<cl_command_queue> gpus_queue;
	std::vector<cl_command_queue> cpus_queue;
	
	std::vector<std::vector<cl_device_id>> gpus;
	std::vector<std::vector<cl_device_id>> cpus;

	std::vector<cl_context> gpus_contexts;
	std::vector<cl_context> cpus_contexts;

	std::vector<cl_platform_id> cpu_platforms;
	std::vector<cl_platform_id> gpu_platforms;

	public:

	Context();
	~Context();

	cl_context getComputeContext();
	cl_command_queue getComputeQueue();

	static Context* initContext();
	static void stopContext();
};

}
}
}
