#ifndef KARU_COMPUTE_ENGINE_H
#define KARU_COMPUTE_ENGINE_H

#include "OpenCL/cl.h"
#include "OpenCL/opencl.h"


class compute_engine {
    cl_context context;
    cl_command_queue queue;

public:
    compute_engine();
    void load_kernel(const char* kernel);
};

#endif