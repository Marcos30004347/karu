#include "algebra/core/compute/ComputeProgram.hpp"
#include "algebra/core/compute/ComputeContext.hpp"

#include <iostream>

using namespace karu;
using namespace algebra;
using namespace compute;

extern Context* karu_core_global_ctx;

Program::Program(const char* src)
{
	cl_int err;
	this->cp_src = src;

	this->cp_program = clCreateProgramWithSource(karu_core_global_ctx->getComputeContext(), 1, (const char **)&src, NULL, &err);
	
	clHandleError(err);

	err = clBuildProgram(this->cp_program, 0, NULL, NULL, NULL, NULL);
	if(err != CL_SUCCESS)
	{
		size_t len;
		err = clGetProgramBuildInfo(this->cp_program, karu_core_global_ctx->getGpuDevice(), CL_PROGRAM_BUILD_LOG, 0, NULL, &len);
		char *buffer = (char*)malloc(len*sizeof(char));
		err = clGetProgramBuildInfo(this->cp_program, karu_core_global_ctx->getGpuDevice(), CL_PROGRAM_BUILD_LOG, len, buffer, NULL);
		std::cout << buffer << std::endl;
		abort();
	}

	clHandleError(err);
}

Program::~Program()
{
	cl_int err = clReleaseProgram(this->cp_program);
	clHandleError(err);
}

cl_program Program::getProgram()
{
	return this->cp_program;
}
