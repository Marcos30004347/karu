#pragma once

#define CL_TARGET_OPENCL_VERSION 220

#ifdef __linux__
#include <CL/cl.h>
#include <CL/opencl.h>
#else
#include "OpenCL/cl.h"
#include "OpenCL/opencl.h"
#endif
// namespace karu {
// namespace algebra {
// namespace compute { 
// namespace cl {


// }
// }
// }
// }
