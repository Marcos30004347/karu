#pragma once

#include <algorithm>
#include <initializer_list>
#include "algebra/core/types.hpp"

namespace karu::algebra
{

class Vector {
	u32  v_size;
	f32* v_data;
public:

	Vector(std::initializer_list<f32> vals);
	Vector(size_t size, f32* data);
	Vector(const Vector& other);
	Vector();
	~Vector();
	const u32 size() const;
	f32* data() const;
	f32& operator[](u32 i);
	
	bool operator ==(Vector& other);
	bool operator ==(Vector&& other);
	Vector& operator =(const Vector& other);

};

}
