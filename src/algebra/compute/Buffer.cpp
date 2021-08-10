#include "algebra/compute/Buffer.hpp"
#include "algebra/compute/Context.hpp"

#include <string.h>
#include <iostream>

using namespace karu;
using namespace algebra;
using namespace compute;

extern Context* karu_core_global_ctx;

Buffer::Buffer(u64 size, type kind, state state)
{
	this->s_kind = kind;
	this->s_size = size;

	if(state == state::MEM_GPU)
	{
		cl_int err;
		switch (this->s_kind)
		{
		case type::READ_WRITE:
			this->s_compute_unit_ref = clCreateBuffer(
				karu_core_global_ctx->getComputeContext(),
				CL_MEM_READ_WRITE,
				this->s_size,
				NULL,
				&err
			);
		
			break;
		case type::READ_ONLY:
			this->s_compute_unit_ref = clCreateBuffer(
				karu_core_global_ctx->getComputeContext(),
				CL_MEM_READ_ONLY,
				this->s_size,
				NULL,
				&err
			);
			break;
		case type::WRITE_ONLY:
			this->s_compute_unit_ref = clCreateBuffer(
				karu_core_global_ctx->getComputeContext(),
				CL_MEM_WRITE_ONLY,
				this->s_size,
				NULL,
				&err
			);
			break;
		default:
			this->s_compute_unit_ref = clCreateBuffer(
				karu_core_global_ctx->getComputeContext(),
				CL_MEM_READ_WRITE,
				this->s_size,
				NULL,
				&err
			);
			break;
		}
		
		clHandleError(err);
	
		this->s_is_compute_unit_allocated = true;
		this->s_is_logic_unit_allocated = false;
	}

	if(state == state::MEM_CPU)
	{
		this->s_logic_unit_ref = malloc(this->s_size);
		this->s_is_logic_unit_allocated = true;
		this->s_should_cleanup = true;
		this->s_is_compute_unit_allocated = false;
	}
}

Buffer::Buffer(void* const data, u64 size, type kind, bool cleanup)
{
		this->s_kind = kind;
		this->s_is_logic_unit_allocated = true;
		this->s_is_compute_unit_allocated = false;
		this->s_size = size;

		this->s_should_cleanup = cleanup;
		this->s_logic_unit_ref = data;
}

Buffer::~Buffer()
{
		if(this->s_should_cleanup && this->s_is_logic_unit_allocated)
		{
			free(this->s_logic_unit_ref);
		}
	
		if(this->s_is_compute_unit_allocated)
		{
			cl_int err;
			err = clReleaseMemObject(this->s_compute_unit_ref);
			clHandleError(err);
		}
}

void* Buffer::download()
{
	this->toLogicUnit();
	return this->s_logic_unit_ref;
}

void Buffer::toComputeUnit()
{
	if(!s_is_compute_unit_allocated)
	{
		cl_int err;
		switch (this->s_kind)
		{
		case type::READ_WRITE:
			this->s_compute_unit_ref = clCreateBuffer(
				karu_core_global_ctx->getComputeContext(),
				CL_MEM_READ_WRITE,
				this->s_size,
				NULL,
				&err
			);
			break;
		case type::READ_ONLY:
			this->s_compute_unit_ref = clCreateBuffer(
				karu_core_global_ctx->getComputeContext(),
				CL_MEM_READ_ONLY,
				this->s_size,
				NULL,
				&err
			);
			break;
		case type::WRITE_ONLY:
			this->s_compute_unit_ref = clCreateBuffer(
				karu_core_global_ctx->getComputeContext(),
				CL_MEM_WRITE_ONLY,
				this->s_size,
				NULL,
				&err
			);
			break;
		default:
			this->s_compute_unit_ref = clCreateBuffer(
				karu_core_global_ctx->getComputeContext(),
				CL_MEM_READ_WRITE,
				this->s_size,
				NULL,
				&err
			);
			break;
		}
	
		clHandleError(err);
	}

	clEnqueueWriteBuffer(
		karu_core_global_ctx->getComputeQueue(),
		this->s_compute_unit_ref,
		CL_TRUE,
		0,
		this->s_size,
		this->s_logic_unit_ref,
		0,
		NULL,
		NULL
	);

	if(this->s_should_cleanup && this->s_is_logic_unit_allocated)
	{
		free(this->s_logic_unit_ref);
	}

	this->s_is_logic_unit_allocated = false;
	this->s_is_compute_unit_allocated = true;
}

void Buffer::toLogicUnit()
{
	if(this->s_is_compute_unit_allocated)
	{
		if(!s_is_logic_unit_allocated)
			this->s_logic_unit_ref = malloc(this->s_size);
			// std::cout << this->s_size << std::endl;
			cl_int err = clEnqueueReadBuffer(
				karu_core_global_ctx->getComputeQueue(),
				this->s_compute_unit_ref,
				CL_TRUE,
				0,
				this->s_size,
				this->s_logic_unit_ref,
				0,
				NULL,
				NULL
			);
		
		clHandleError(err);

		err = clReleaseMemObject(this->s_compute_unit_ref);
		
		clHandleError(err);
	
		this->s_is_logic_unit_allocated = true;
		this->s_is_compute_unit_allocated = false;
	}
}

void* Buffer::upload()
{
	if(!this->s_is_compute_unit_allocated)
	{
		this->toComputeUnit();
	}

	return (void*)&this->s_compute_unit_ref;
}
