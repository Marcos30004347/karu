#pragma once

#include "algebra/core/types.hpp"
#include "algebra/core/compute/OpenCL.hpp"

namespace karu { 
namespace algebra { 
namespace compute {

class Storage {
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

	Storage(i8* data, u64 size, type kind = type::READ_WRITE, bool copy = false);
	Storage(u64 size, type kind, state state);

	~Storage();

	void* data();

	void* ref();

	void toComputeUnit(bool free_cpu_data = false);
	void toLogicUnit(bool free_gpu_data = false);

	private:


	bool s_is_logic_unit_allocated;
	bool s_is_compute_unit_allocated;

	state   		s_state;
	cl_mem  s_compute_unit_ref;
	i8*     		s_logic_unit_ref;
	u64     		s_size;
	type    		s_kind;

};

}
}
}
