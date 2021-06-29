#include "compute.hpp"

#include "OpenCL/cl.h"
#include "OpenCL/opencl.h"

#include <stdio.h>
#include <fstream>
#include <string>
#include <time.h>

#include "types.hpp"

compute_engine::compute_engine() {

	std::ifstream ifs("src/kernels/matrix_mul.cl");
  	std::string content((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()) );

	cl_int err;

	this->context = clCreateContextFromType(NULL, CL_DEVICE_TYPE_GPU, NULL, NULL, &err);
    
    cl_uint num = 1;

	clGetDeviceIDs(0, CL_DEVICE_TYPE_GPU, 0, NULL, &num);

	cl_device_id devices[num];

	clGetDeviceIDs(0, CL_DEVICE_TYPE_GPU, num, devices, &num);

    this->queue = clCreateCommandQueue(context, devices[0], 0, &err);
};

void compute_engine::load_kernel(const char* kernel) {
    std::ifstream ifs("src/kernels/matrix_mul.cl");
  	std::string content((std::istreambuf_iterator<char>(ifs)),(std::istreambuf_iterator<char>()    ) );
	
    const char* source = content.c_str();

}



