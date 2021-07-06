#include "algebra/core/compute/ComputeProgram.hpp"
#include "algebra/core/compute/ComputeContext.hpp"

using namespace karu;
using namespace algebra;
using namespace compute;

extern Context* karu_core_global_ctx;

Program::Program(const char* src)
{
	cl_int err;
	this->cp_src = src;
	this->cp_program = clCreateProgramWithSource(karu_core_global_ctx->getComputeContext(), 1, (const char **)&src, NULL, &err);
	err = clBuildProgram(this->cp_program, 0, NULL, NULL, NULL, NULL);
}

Program::~Program()
{
	clReleaseProgram(this->cp_program);
}

cl_program Program::getProgram()
{
	return this->cp_program;
}
