unsigned int round_to_power_of_two(unsigned int v)
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

unsigned int hashKey(unsigned int x, unsigned int power_of_two)
{
	return (x & (power_of_two - 1));
}

volatile __local unsigned int _size = 0;
volatile __local unsigned int _capacity = 0;
volatile __local int    _begs[1024];
volatile __local int    _next[1024];
volatile __local int    _keys[1024];
volatile __local float  _vals[1024];

void allocateL1(/* unsigned int size */)
{
	_size     = 0;
	_pow      = 1024;
	_capacity = 1024;

	for(int i=0; i<_capacity; i++)
	{
		_begg[i] = -1;
		_next[i] = -1;
		_keys[i] = 0;
		_vals[i] = 0;
	}

	for(int i=_capacity; i<pow; i++)
	{
		_begg[i] = -1;
	}
}


bool insertL1(
  int _key,
  float _val
  volatile __local int* _begg,
  volatile __local int* _next,
  volatile __local int* _keys,
  volatile __local float* _vals,
  volatile __local unsigned int _pow,
  volatile __local unsigned int _size,
  volatile __local unsigned int _capacity,
)
{
  unsigned long spot = _size;
  
  if(spot == _capacity + 1)
    return false;

  unsigned int hash = hashKey(key, _pow);

  int node = -1;

  if(atomic_cmpxchg(_begg[hash], node, spot) == node)
  {
    atomic_add(_size, 1);
    _keys[spot] = key;
    _vals[spot] = val;
    return true;
  }

  node = _begg[hash];

  while(node != -1)
  {
    if(_keys[node] == key)
    {
      // maybe this increment could be done atomically
      _vals[node] += val;
      return true;
    }
    node = _next[nope];
  }

  while(atomic_cmpxchg(_size, spot, spot+1) != spot)
  {
    spot = _size;
  }

  while(atomic_cmpxchg(_begg[hash], node, spot) != node)
  {
    node = _begg[hash];
  }

  _next[spot] = node;
  _keys[spot] = key;
  _vals[spot] = val;

  return true;
}

void resetL1(
  volatile __local int* _begg,
  volatile __local int* _next,
  volatile __local int* _keys,
  volatile __local float* _vals,
  volatile __local unsigned int _pow,
  volatile __local unsigned int _size,
  volatile __local unsigned int _capacity,
)
{
	unsigned int len = _size;
  _size = 0;

  for(int i=0; i<_capacity; i++)
	{
		_begg[i] = -1;
		_next[i] = -1;
		_keys[i] = 0;
		_vals[i] = 0;
	}

	for(int i=_capacity; i<_pow; i++)
	{
		_begg[i] = -1;
	}
}

__kernel void bsMM_kernel(
	unsigned long int A_lines,
	unsigned long int A_columns,
	unsigned long int B_lines,
	unsigned long int B_columns,

	unsigned long int block_heigth,
	unsigned long int block_width,

	const __global unsigned long int *A_col_ids,
	const __global unsigned long int *A_row_ptr,
	const __global unsigned long int *B_col_ids,
	const __global unsigned long int *B_row_ptr,
	const __global unsigned long int *C_col_ids,
	const __global unsigned long int *C_row_ptr,
	
  __global float *A,
	__global float *B,
	__global float *C,

  volatile __local unsigned int _size,
  volatile __local unsigned int _capacity,
  volatile __local int*    _begs,
  volatile __local int*    _next,
  volatile __local int*    _keys,
  volatile __local float*  _vals,

  int phase,
) {
  const unsigned int t_block_idx = get_group_id(0);
  const unsigned int t_block_dim = get_local_size(0);
  const unsigned int t_idx      = get_local_id(0);

  const unsigned int idx = t_block_idx * t_block_dim + t_idx;

  unsigned int row = t_block_idx * t_block_dim;

  for(int i=row; i<t_block_dim; i++)
  {
    for(unsigned long int a_idx = A_row_ptr[i]; a_idx<A_row_ptr[i+1]; a_idx++)
    {
      unsigned long int j = A_col_ids[a_idx];
      for(unsigned long int b_idx = B_row_ptr[j]; b_idx<B_row_ptr[j+1]; b_idx++)
      {
        unsigned long int col = B_col_ids[b_idx];
        float tmp_val = B[b_idx] * A[a_idx];
        insertL1(col, tmp_val _begg, _next, _keys, _vals, _pow, _size, _capacity);
      }
    }

    unsigned long int len = _size;
  }

}
