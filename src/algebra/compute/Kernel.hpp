#pragma once

#include <vector>
#include <iostream>

#include "algebra/core/types.hpp"
#include "algebra/compute/OpenCL.hpp"
#include "algebra/compute/Buffer.hpp"
#include "algebra/compute/Program.hpp"

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

	// template<typename T>
	// void setKernelArgument(u32 id, T&& value);

	template<typename T>
	void setKernelArgument(u32 id, T& value)
	{
		cl_int err = clSetKernelArg(this->ck_kernel, id, sizeof(T), &value);
		clHandleError(err);
	}

	template<typename T>
	void setKernelArgument(u32 id, T* value)
	{
		cl_int err = clSetKernelArg(this->ck_kernel, id, sizeof(T), value);
		clHandleError(err);
	}

	template<typename T>
	void setKernelArgument(u32 id, const T* const value)
	{
		cl_int err = clSetKernelArg(this->ck_kernel, id, sizeof(T), value);
		clHandleError(err);
	}

	template<typename T>
	void setKernelArgument(u32 id, const T* value)
	{
		cl_int err = clSetKernelArg(this->ck_kernel, id, sizeof(T), value);
		clHandleError(err);
	}

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
