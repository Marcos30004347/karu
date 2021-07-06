#include "algebra/core/compute/ComputeContext.hpp"

using namespace karu;
using namespace algebra;
using namespace compute;

Context* karu_core_global_ctx = nullptr;

Context* Context::initContext()
{
	karu_core_global_ctx = new Context();
	return karu_core_global_ctx;
}

Context::Context()
{
	cl_int err;

	this->context = clCreateContextFromType(nullptr, CL_DEVICE_TYPE_ALL, nullptr, nullptr, &err);

	cl_uint num = 1;
	cl_device_id devices;

	clGetDeviceIDs(0, CL_DEVICE_TYPE_GPU, 0, nullptr, &num);
	clGetDeviceIDs(0, CL_DEVICE_TYPE_GPU, num, &devices, &num);
	
	this->queue = clCreateCommandQueueWithProperties(context, devices, nullptr, &err);
}

Context::~Context()
{
	clReleaseContext(this->context);
	clReleaseCommandQueue(this->queue);
}

cl_command_queue Context::getComputeQueue()
{
	return this->queue;
}

cl_context Context::getComputeContext()
{
	return this->context;
}


