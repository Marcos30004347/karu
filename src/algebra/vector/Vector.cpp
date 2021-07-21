#include "algebra/vector/Vector.hpp"

namespace karu::algebra
{

Vector::Vector(std::initializer_list<f32> vals)
{
	v_size = vals.size();
	this->v_data = new f32[vals.size()];
	std::copy(vals.begin(), vals.end(), v_data);
}

Vector::Vector(size_t size, f32* data)
{
	v_size = size;
	this->v_data = new f32[size];
	std::copy(&data[0], &data[size-1], v_data);
}

Vector::~Vector()
{
	delete this->v_data;
}

f32* Vector::data()
{
	return this->v_data;
}

f32& Vector::operator[](u32 i)
{
	return this->v_data[i];
}

u32 Vector::size()
{
	return this->v_size;
}

bool Vector::operator ==(Vector& other)
{
	for(i32 i=0; i<this->v_size; i++)
	{
		if(this->v_data[i] != other[i])
			return false;
	}
	return true;
}

bool Vector::operator ==(Vector&& other)
{
	for(i32 i=0; i<this->v_size; i++)
	{
		if(this->v_data[i] != other[i])
			return false;
	}
	return true;
}

}
