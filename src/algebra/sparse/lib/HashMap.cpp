#include "HashMap.hpp"

#include <pthread.h>

using namespace karu;

u32 round_to_power_of_two(u32 v)
{
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;
    return v;
}

HashMap::HashMap(u32 size)
{
	u32 pow = round_to_power_of_two(size);

	this->_size.store(0);
	this->_begg = new std::atomic<i32>[pow];

	this->_next = new i32[size+1];
	this->_keys = new i32[size+1];
	this->_vals = new f32[size+1];

	for(int i=0; i<size; i++)
	{
		this->_begg[i].store(-1);
		this->_next[i] = -1;
		this->_keys[i] = 0;
		this->_vals[i] = 0;
	}

	for(int i=size; i<pow; i++)
	{
		this->_begg[i] = -1;
	}

	this->_pow  = pow;
	this->_capacity = size;
}

HashMap::~HashMap()
{
	delete this->_begg;
	delete this->_next;
	delete this->_keys;
	delete this->_vals;
}

u32 HashMap::hashKey(u32 x, u32 power_of_two)
{
	return (x & (power_of_two - 1));
}

void HashMap::reset()
{
	u32 len = this->_size.load();

	this->_size.store(0);
	
	for(int i=0; i<this->_capacity; i++)
	{
		this->_begg[i].store(-1);
		this->_next[i] = -1;
		this->_keys[i] = 0;
		this->_vals[i] = 0;
	}

	for(int i=this->_capacity; i<this->_pow; i++)
	{
		this->_begg[i] = -1;
	}
}

bool HashMap::insert(i32 key, f32 val)
{
	u32 spot 	= this->_size.load();

	if(spot == this->_capacity + 1)
		return false;

	u32 hash  = hashKey(key, this->_pow);

	i32 node = -1;

	if(this->_begg[hash].compare_exchange_weak(node, spot))
	{
		// empy hash, just insert the node
		this->_size.fetch_add(1);
		this->_keys[spot] = key;
		this->_vals[spot] = val;

		return true;
	}

	// Iterate Linked List nodes to accumulate in existing keys
	node = this->_begg[hash].load();
	while(node != -1)
	{
		if(this->_keys[node] == key)
		{

			this->_vals[node] += val;
			return true;
		}
	
		node = this->_next[node];
	}

	// If no node found in Linked List, add a new one as head
	node = this->_begg[hash].load();
	spot = this->_size.load();

	while(!this->_size.compare_exchange_weak(spot, spot+1))
		spot = this->_size.load();
	
	while(!this->_begg[hash].compare_exchange_weak(node, spot))
	{
		node = this->_begg[hash].load();
	}

	this->_next[spot] = node;
	this->_keys[spot] = key;
	this->_vals[spot] = val;

	return true;
}

f32 HashMap::get(u32 key)
{
	u64 node = this->_begg[this->hashKey(key, this->_pow)];
	
	if(node == -1)
		return 0.0;

	while(key != this->_keys[node])
	{
		node = this->_next[node];
		
		if(node == -1)
			return 0.0;
	}

	return this->_vals[node];
}

u32 HashMap::size()
{
	return this->_size.load();
}

u32 HashMap::capacity() {
	return this->_capacity;
}
