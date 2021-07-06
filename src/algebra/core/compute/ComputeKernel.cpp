#include <assert.h>

#include "algebra/core/compute/ComputeKernel.hpp"
#include "algebra/core/compute/ComputeContext.hpp"

using namespace karu;
using namespace algebra;
using namespace compute;

extern Context* karu_core_global_ctx;

Event::Event(){}

Event::~Event(){
	clReleaseEvent(this->e_event);
}

cl_event* Event::ref()
{
	return &this->e_event;
}
Kernel::Kernel(Program* program, const char* kernel_name)
{
	this->ck_kernel = clCreateKernel(program->getProgram(), kernel_name, NULL);
}

Kernel::~Kernel()
{
	clReleaseKernel(this->ck_kernel);
}

void Kernel::setKernelArgument(u32 id, u32 size, void* ptr)
{
	clSetKernelArg(this->ck_kernel, id, size, ptr);
}

void Kernel::enqueue(std::vector<u64> global_work_size, std::vector<u64> local_work_size, std::vector<Event> wait_list, Event* event)
{
	assert(global_work_size.size() == local_work_size.size());

	cl_int err;
	err = clEnqueueNDRangeKernel(karu_core_global_ctx->getComputeQueue(), this->ck_kernel, global_work_size.size(), NULL, global_work_size.data(), local_work_size.data(), wait_list.size(), (cl_event*)wait_list.data(), event->ref());
}
