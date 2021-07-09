#include "algebra/core/compute/ComputeContext.hpp"
#include "algebra/matrix/kernels/kernels.hpp"


#include <iostream>
#include <string.h>
using namespace karu;
using namespace algebra;
using namespace compute;

Context* karu_core_global_ctx = nullptr;

Context* Context::initContext()
{
	karu_core_global_ctx = new Context();
	create_sparse_matrix_kernels();
	return karu_core_global_ctx;
}

void Context::stopContext()
{
	destroy_sparse_matrix_kernels();
	delete karu_core_global_ctx;
}

Context::Context()
{
	cl_int err;
	cl_uint numGPUs = 0;
	cl_uint numCPUs = 0;
	cl_uint numPlatforms = 1;

	err = clGetPlatformIDs(0, NULL, &numPlatforms);
	clHandleError(err);

	cl_platform_id platforms[numPlatforms];
	err = clGetPlatformIDs(numPlatforms, platforms, NULL);
	clHandleError(err);

	for(int i=0; i<numPlatforms; i++)
	{
		numGPUs = 0;
		numCPUs = 0;
	
		clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_GPU, 0, nullptr, &numGPUs);
	
		if(numGPUs != 0)
		{
			std::vector<cl_device_id> gpus = std::vector<cl_device_id>(numGPUs, nullptr);

			this->gpu_platforms.push_back(platforms[i]);

			err = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_GPU, numGPUs, gpus.data(), &numGPUs);
			clHandleError(err);

			this->gpus.push_back(gpus);
		}
	
		clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_CPU, 0, nullptr, &numCPUs);
		if(numCPUs != 0)
		{
			std::vector<cl_device_id> cpus = std::vector<cl_device_id>(numCPUs, nullptr);
		
			this->cpu_platforms.push_back(platforms[i]);
			
			err = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_CPU, numCPUs, cpus.data(), &numCPUs);
			clHandleError(err);
			
			this->cpus.push_back(cpus);
		}
	}

	
	for(int i=0; i<this->gpu_platforms.size(); i++)
	{
		for(cl_device_id gpu: this->gpus[i])
		{
			cl_context_properties properties[] = { CL_CONTEXT_PLATFORM, (cl_context_properties)gpu_platforms[i], 0 };
			cl_context context = clCreateContextFromType(properties, CL_DEVICE_TYPE_ALL, NULL, NULL, &err);
			clHandleError(err);
		
			this->gpus_contexts.push_back(context);
		
			cl_command_queue queue = clCreateCommandQueueWithProperties(context, gpu, nullptr, &err);
			clHandleError(err);

			this->gpus_queue.push_back(queue);
		}
	}

	for(int i=0; i<this->cpu_platforms.size(); i++)
	{
		for(cl_device_id cpu: this->cpus[i])
		{
			cl_context_properties properties[] = { CL_CONTEXT_PLATFORM, (cl_context_properties)cpu_platforms[i], 0 };
			cl_context context = clCreateContextFromType(properties, CL_DEVICE_TYPE_ALL, NULL, NULL, &err);
			clHandleError(err);
		
			this->cpus_contexts.push_back(context);
		
			cl_command_queue queue = clCreateCommandQueueWithProperties(context, cpu, nullptr, &err);
			clHandleError(err);

			this->cpus_queue.push_back(queue);
		}
	}
}

Context::~Context()
{
	cl_int err;

	
	for(int i=0; i<this->gpu_platforms.size(); i++)
	{
		for(cl_device_id gpu: this->gpus[i])
		{
			clReleaseDevice(gpu);
		}
	}

	for(int i=0; i<this->cpu_platforms.size(); i++)
	{
		for(cl_device_id cpu: this->cpus[i])
		{
			clReleaseDevice(cpu);
		}
	}

	for(cl_context c : this->gpus_contexts)
	{
		err = clReleaseContext(c);
		clHandleError(err);
	}

	for(cl_context c : this->cpus_contexts)
	{
		err = clReleaseContext(c);
		clHandleError(err);
	}

	for(cl_command_queue q : this->gpus_queue)
	{
		err = clReleaseCommandQueue(q);
		clHandleError(err);
	}

	for(cl_command_queue q : this->cpus_queue)
	{
		err = clReleaseCommandQueue(q);
		clHandleError(err);
	}
}

cl_command_queue Context::getComputeQueue()
{
	return this->gpus_queue[0];
}

cl_context Context::getComputeContext()
{
	return this->gpus_contexts[0];
}


