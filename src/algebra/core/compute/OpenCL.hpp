#pragma once

#define CL_TARGET_OPENCL_VERSION 220

#ifdef __linux__
#include <CL/cl.h>
#include <CL/opencl.h>
#else
#include "OpenCL/cl.h"
#include "OpenCL/opencl.h"
#endif

struct clPlatformInfo 
{
	char name[24];
	char vendor[24];
	char version[24];
	char profile[24];
	char extensions[1024];
};

const char *getErrorString(cl_int error);

void clHandleError(cl_int result);

clPlatformInfo getPlatformInfo(cl_platform_id id);
