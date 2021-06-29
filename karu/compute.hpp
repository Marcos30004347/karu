#ifndef KARU_COMPUTE_ENGINE_H
#define KARU_COMPUTE_ENGINE_H

#ifdef __linux__
#include <CL/cl.h>
#include <CL/opencl.h>
#else
#include "OpenCL/cl.h"
#include "OpenCL/opencl.h"
#endif

#include <string.h>

#include "types.hpp"
namespace karu {

// TODO: initialize queue
cl_command_queue queue;

// TODO: initialize context
cl_context context;

class Memory {
	enum type;

	enum state {
		MEM_CPU = 0,
		MEM_GPU,
	};

	bool _is_cpu_allocated;
	bool _is_gpu_allocated;

	state   _state;
	cl_mem  _gpu_ref;
	i8*     _cpu_ref;
	u64     _size;
	type    _kind;

public:
	enum type {
		READ_ONLY = 0,
		WRITE_ONLY,
		READ_WRITE,
	};

	Memory(i8* data, u64 size, type kind = type::READ_WRITE, bool copy = true)
	{
		this->_kind = kind;
		this->_is_cpu_allocated = true;
		this->_is_gpu_allocated = false;
	
		this->_state = state::MEM_CPU;
		if(!copy)
		{
			this->_size = size;
			this->_cpu_ref = data;
		}
		else
		{
			this->_cpu_ref = (i8*)malloc(size);
			this->_size = size;
			memcpy(this->_cpu_ref, data, size);
		}
	}

	~Memory()
	{
		if(this->_is_cpu_allocated)
			free(this->_cpu_ref);
		if(this->_is_gpu_allocated)
			clReleaseMemObject(this->_gpu_ref);
	}

	void* data()
	{
		if(this->_state == state::MEM_GPU)
		{
			clEnqueueReadBuffer(
				queue,
				this->_gpu_ref, 
				CL_TRUE,
				0,
				this->_size,
				(void*)this->_cpu_ref,
				0,
				NULL,
				NULL
			);
		}
	}

	void toGPU(bool free_cpu_data = false)
	{
		if(!_is_gpu_allocated)
		{
			switch (this->_kind)
			{
			case type::READ_WRITE:
				this->_gpu_ref = clCreateBuffer(context, CL_MEM_READ_WRITE, this->_size, NULL, NULL),
				break;
			case type::READ_ONLY:
				this->_gpu_ref = clCreateBuffer(context, CL_MEM_READ_ONLY, this->_size, NULL, NULL),
			case type::WRITE_ONLY:
				this->_gpu_ref = clCreateBuffer(context, CL_MEM_WRITE_ONLY, this->_size, NULL, NULL),
			default:
				this->_gpu_ref = clCreateBuffer(context, CL_MEM_READ_WRITE, this->_size, NULL, NULL),
				break;
			}
			
			this->_is_gpu_allocated = true;
		}
	
		clEnqueueWriteBuffer(queue, this->_gpu_ref, CL_TRUE, 0, this->_size, this->_cpu_ref, 0, NULL, NULL);

		if(free_cpu_data)
		{
			free(this->_cpu_ref);
			this->_is_cpu_allocated = false;
		}
	}

	void toCPU(bool free_gpu_data = false)
	{
		if(!_is_cpu_allocated)
		{
			this->_cpu_ref = (i8*)malloc(this->_size);
			this->_is_cpu_allocated = true;
		}
	
		// TODO: handle err
		cl_int err = clEnqueueReadBuffer(queue, this->_gpu_ref, CL_TRUE, 0, this->_size, this->_cpu_ref, 0, NULL, NULL);

		if(free_gpu_data)
		{
			clReleaseMemObject(this->_gpu_ref);
			this->_is_gpu_allocated = false;
		}
	}
}

}

class compute_engine {
    cl_context context;
    cl_command_queue queue;
};

#endif
