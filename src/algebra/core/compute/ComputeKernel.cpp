#include <assert.h>

#include "algebra/core/compute/ComputeKernel.hpp"
#include "algebra/core/compute/ComputeContext.hpp"

using namespace karu;
using namespace algebra;
using namespace compute;

extern Context* karu_core_global_ctx;

Event::Event(){}

Event::~Event(){
	cl_int err = clReleaseEvent(this->e_event);
	clHandleError(err);
}

cl_event* Event::ref()
{
	return &this->e_event;
}
Kernel::Kernel(Program* program, const char* kernel_name)
{
	cl_int err;
	this->ck_kernel = clCreateKernel(program->getProgram(), kernel_name, &err);
	clHandleError(err);
}

Kernel::~Kernel()
{
	cl_int err = clReleaseKernel(this->ck_kernel);
	clHandleError(err);
}

void Kernel::setKernelArgument(u32 id, u32 size, void* ptr)
{
	cl_int err = clSetKernelArg(this->ck_kernel, id, size, ptr);
	clHandleError(err);
}

void Kernel::enqueue(std::vector<u64> global_work_size, std::vector<u64> local_work_size, std::vector<Event> wait_list, Event* event)
{
	assert(global_work_size.size() == local_work_size.size());

	cl_int err;
	err = clEnqueueNDRangeKernel(karu_core_global_ctx->getComputeQueue(), this->ck_kernel, global_work_size.size(), NULL, global_work_size.data(), local_work_size.data(), wait_list.size(), (cl_event*)wait_list.data(), event->ref());
	clHandleError(err);
}

void Kernel::enqueue(std::vector<u64> global_work_size, std::vector<u64> local_work_size)
{
	assert(global_work_size.size() == local_work_size.size());

	cl_int err;
	err = clEnqueueNDRangeKernel(karu_core_global_ctx->getComputeQueue(), this->ck_kernel, global_work_size.size(), NULL, global_work_size.data(), local_work_size.data(), 0, nullptr, nullptr);
	clHandleError(err);
}
