#pragma once

#include <atomic>
#include "algebra/core/types.hpp"

namespace karu {

class HashMap
{
public:
	bool insert(i32 key, f32 val);
	void reset();

	f32 get(u32 key);

	u32 size();
	u32 capacity();

	HashMap(u32 capacity);
	~HashMap();

	i32* _keys;
	f32* _vals;
	i32 	_pow;
	i32* _next;

	std::atomic<i32>* _begg;
	std::atomic<u32> _size;

	i32  _capacity;
private:
	u32 hashKey(u32 x, u32 power_of_two);
};

}

