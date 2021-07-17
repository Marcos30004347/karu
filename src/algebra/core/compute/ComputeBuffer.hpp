#pragma once

#include "algebra/core/types.hpp"
#include "algebra/core/compute/OpenCL.hpp"

namespace karu { 
namespace algebra { 
namespace compute {

#define BUFFER_ARG_SIZE sizeof(cl_mem)

class Buffer {
	friend class Kernel;

	public:
	enum type {
		READ_ONLY = 0,
		WRITE_ONLY,
		READ_WRITE,
	};

	enum state {
		MEM_CPU = 0,
		MEM_GPU,
	};

	Buffer(void* buffer, u64 size, type kind, bool cleanup);
	Buffer(u64 size, type kind, state state);

	~Buffer();

	void* download();
	void* upload();
	
	void toComputeUnit();
	void toLogicUnit();

	private:
	bool s_is_logic_unit_allocated;
	bool s_is_compute_unit_allocated;
	
	bool s_should_cleanup;
	
	cl_mem  s_compute_unit_ref;
	void*     		s_logic_unit_ref;
	u64     		s_size;
	type    		s_kind;
};

}
}
}
