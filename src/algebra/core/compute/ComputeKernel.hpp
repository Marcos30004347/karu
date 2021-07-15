#pragma once

#include <vector>

#include "algebra/core/types.hpp"
#include "algebra/core/compute/OpenCL.hpp"
#include "algebra/core/compute/ComputeBuffer.hpp"
#include "algebra/core/compute/ComputeProgram.hpp"

namespace karu {
namespace algebra {
namespace compute {

class Event {
	public:
	Event();
	~Event();

	cl_event* ref();

	private:
	cl_event e_event;
};

class Kernel {
	public:

	Kernel(Program* program, const char* kernel_name);
	~Kernel();

	void setKernelArgument(u32 id, u32 size, void* ptr);
	
	void enqueue(std::vector<u64> global_work_size, std::vector<u64> local_work_size, std::vector<Event> wait_list, Event* event);
	void enqueue(std::vector<u64> global_work_size, std::vector<u64> local_work_size);

	size_t getWorkGroupSize();

	private:
	cl_kernel ck_kernel;
};

}
}
}
