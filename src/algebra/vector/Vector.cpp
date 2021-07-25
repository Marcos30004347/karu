#include "algebra/vector/Vector.hpp"
#include <iostream>
namespace karu::algebra
{
Vector::Vector()
{
	v_size = 0;
	this->v_data = nullptr;
}

Vector::Vector(std::initializer_list<f32> vals)
{
	v_size = vals.size();
	this->v_data = new f32[vals.size()];
	std::copy(vals.begin(), vals.end(), v_data);
}


Vector::Vector(const Vector& other)
{
	v_size = other.v_size;
	this->v_data = new f32[other.v_size];
	std::copy(other.v_data, other.v_data + other.v_size, v_data);
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

f32* Vector::data() const
{
	return this->v_data;
}

f32& Vector::operator[](u32 i)
{
	return this->v_data[i];
}

const u32 Vector::size() const
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

Vector& Vector::operator =(const Vector& other)
{
	v_size = other.v_size;
	this->v_data = new f32[other.v_size];
	std::copy(other.v_data, other.v_data + other.v_size, v_data);
	return *this;
}


}
