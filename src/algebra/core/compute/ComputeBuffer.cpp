#include "algebra/core/compute/ComputeBuffer.hpp"
#include "algebra/core/compute/ComputeContext.hpp"

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
	this->s_state = state;

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
		this->s_logic_unit_ref = (i8*)malloc(this->s_size);
		this->s_is_logic_unit_allocated = true;
		this->s_is_compute_unit_allocated = false;
	}
}

Buffer::Buffer(void* data, u64 size, type kind, bool copy)
{
		this->s_kind = kind;
		this->s_is_logic_unit_allocated = true;
		this->s_is_compute_unit_allocated = false;
		this->s_size = size;
	
		this->s_state = state::MEM_CPU;
	
		if(!copy)
		{
			this->s_logic_unit_ref = data;
		}
		else
		{
			this->s_logic_unit_ref = (i8*)malloc(size);
			memcpy(this->s_logic_unit_ref, data, size);
		}
}

Buffer::~Buffer()
{
		if(this->s_is_logic_unit_allocated)
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

void* Buffer::data()
{
		if(this->s_state == state::MEM_GPU)
		{
			cl_int err = clEnqueueReadBuffer(
				karu_core_global_ctx->getComputeQueue(),
				this->s_compute_unit_ref, 
				CL_TRUE,
				0,
				this->s_size,
				(void*)this->s_logic_unit_ref,
				0,
				NULL,
				NULL
			);

			clHandleError(err);
		}
	
		return this->s_logic_unit_ref;
}

void Buffer::toComputeUnit(bool free_cpu_data)
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
	
		this->s_is_compute_unit_allocated = true;
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

	if(free_cpu_data)
	{
		free(this->s_logic_unit_ref);
		this->s_is_logic_unit_allocated = false;
	}
}

void Buffer::toLogicUnit(bool free_gpu_data)
{
	if(!s_is_logic_unit_allocated)
	{
		this->s_logic_unit_ref = (i8*)malloc(this->s_size);
		this->s_is_logic_unit_allocated = true;
	}

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

	if(free_gpu_data)
	{
		cl_int err = clReleaseMemObject(this->s_compute_unit_ref);
		
		clHandleError(err);
		
		this->s_is_compute_unit_allocated = false;
	}
}

void* Buffer::ref()
{
	if(this->s_is_compute_unit_allocated)
	{
		return (void*)&this->s_compute_unit_ref;
	}
	else if(this->s_is_logic_unit_allocated)
	{
		return (void*)&this->s_logic_unit_ref;
	}

	return nullptr;
}
