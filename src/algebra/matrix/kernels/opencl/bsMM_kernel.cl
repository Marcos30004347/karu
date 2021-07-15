// unsigned int round_to_power_of_two(
//   unsigned int v
// )
// {
//     v--;
//     v |= v >> 1;
//     v |= v >> 2;
//     v |= v >> 4;
//     v |= v >> 8;
//     v |= v >> 16;
//     v++;
//     return v;
// }

// unsigned int hashKey(
//   unsigned int x,
//   unsigned int power_of_two
// )
// {
// 	return (x & (power_of_two - 1));
// }

// void allocateL1(
//   volatile __local int* _begg,
//   volatile __local int* _next,
//   volatile __local int* _keys,
//   volatile __local float* _vals,
//   volatile __local unsigned int _pow,
//   volatile __local unsigned int _size,
//   volatile __local unsigned int _capacity
// )
// {
// 	_size     = 0;
// 	_pow      = 1024;
// 	_capacity = 1024;

// 	for(int i=0; i<_capacity; i++)
// 	{
// 		_begg[i] = -1;
// 		_next[i] = -1;
// 		_keys[i] = 0;
// 		_vals[i] = 0;
// 	}

// 	for(int i=_capacity; i<pow; i++)
// 	{
// 		_begg[i] = -1;
// 	}
// }


// bool insertL1(
//   int _key,
//   float _val
//   volatile __local int* _begg,
//   volatile __local int* _next,
//   volatile __local int* _keys,
//   volatile __local float* _vals,
//   volatile __local unsigned int _pow,
//   volatile __local unsigned int _size,
//   volatile __local unsigned int _capacity
// )
// {
//   unsigned long spot = _size;
  
//   if(spot == _capacity + 1)
//     return false;

//   unsigned int hash = hashKey(key, _pow);

//   int node = -1;

//   if(atomic_cmpxchg(_begg[hash], node, spot) == node)
//   {
//     atomic_add(_size, 1);
//     _keys[spot] = key;
//     _vals[spot] = val;
//     return true;
//   }

//   node = _begg[hash];

//   while(node != -1)
//   {
//     if(_keys[node] == key)
//     {
//       // maybe this increment could be done atomically
//       _vals[node] += val;
//       return true;
//     }
//     node = _next[nope];
//   }

//   while(atomic_cmpxchg(_size, spot, spot+1) != spot)
//   {
//     spot = _size;
//   }

//   while(atomic_cmpxchg(_begg[hash], node, spot) != node)
//   {
//     node = _begg[hash];
//   }

//   _next[spot] = node;
//   _keys[spot] = key;
//   _vals[spot] = val;

//   return true;
// }

// void resetL1(
//   volatile __local int* _begg,
//   volatile __local int* _next,
//   volatile __local int* _keys,
//   volatile __local float* _vals,
//   volatile __local unsigned int _pow,
//   volatile __local unsigned int _size,
//   volatile __local unsigned int _capacity
// )
// {
// 	unsigned int len = _size;
//   _size = 0;

//   for(int i=0; i<_capacity; i++)
// 	{
// 		_begg[i] = -1;
// 		_next[i] = -1;
// 		_keys[i] = 0;
// 		_vals[i] = 0;
// 	}

// 	for(int i=_capacity; i<_pow; i++)
// 	{
// 		_begg[i] = -1;
// 	}
// }


// unsigned int getMax(
//   volatile __local unsigned int* keys,
//   unsigned int n
// )
// {
//     unsigned int mx = keys[0];
//     for (u32 i = 1; i < n; i++)
//         if (keys[i] > mx)
//             mx = keys[i];
//     return mx;
// }

// void countSort(
//   volatile __local unsigned int* keys,
//   volatile __local float* vals,
//   volatile __local unsigned int* out_keys,
//   volatile __local float* out_vals,
//   unsigned int n,
//   unsigned int exp
// )
// {
//     int i; 
//     int count[10];
    
//     for (i = 1; i < 10; i++)
//       count[i] = 0;
    
//     for (i = 0; i < n; i++)
//         count[(keys[i] / exp) % 10]++;
 
  
//     for (i = 1; i < 10; i++)
//         count[i] += count[i - 1];

//     for (i = n - 1; i >= 0; i--) {
//         out_keys[count[(keys[i] / exp) % 10] - 1] = keys[i];
//         out_vals[count[(keys[i] / exp) % 10] - 1] = vals[i];
//         count[(keys[i] / exp) % 10]--;
//     }
 
//     for (i = 0; i < n; i++)
// 		{
// 			keys[i] = out_keys[i];
// 			vals[i] = out_vals[i];
// 		}
// }

// void sortKeyValuePair(
//   volatile __local unsigned int* keys,
//   volatile __local float* vals,
//   volatile __local unsigned int* out_keys,
//   volatile __local float* out_vals,
//   unsigned int n
// )
// {
//   unsigned int m = getMax(keys, n);
//   for (int exp = 1; m / exp > 0; exp *= 10)
//     countSort(keys, vals, out_keys, out_vals, n, exp);
// }

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

  volatile __local unsigned int size,
  volatile __local unsigned int capacity,
  volatile __local int*    begs,
  volatile __local int*    next,
  volatile __local int*    keys,
  volatile __local float*  vals,
  volatile __local int*    out_keys,
  volatile __local float* out_vals,  
  int phase
) {
  // const unsigned int t_block_idx = get_group_id(0);
  // const unsigned int t_block_dim = get_local_size(0);
  // const unsigned int t_idx      = get_local_id(0);

  // const unsigned int idx = t_block_idx * t_block_dim + t_idx;

  // unsigned int row = t_block_idx * t_block_dim;

  // allocateL1(_begg, _next, _keys, _vals, _pow, _size, _capacity);

  // for(int i=row; i<t_block_dim; i++)
  // {
  //   for(unsigned long int a_idx = A_row_ptr[i]; a_idx<A_row_ptr[i+1]; a_idx++)
  //   {
  //     unsigned long int j = A_col_ids[a_idx];
  //     for(unsigned long int b_idx = B_row_ptr[j]; b_idx<B_row_ptr[j+1]; b_idx++)
  //     {
  //       unsigned long int col = B_col_ids[b_idx];
  //       float tmp_val = B[b_idx] * A[a_idx];
  //       insertL1(col, tmp_val _begg, _next, _keys, _vals, _pow, _size, _capacity);
  //     }
  //   }

  //   unsigned long int len = _size;
  
  //   if(phase == 0)
  //   {
  //     C_row_ptr[i] = len;
  //   }
  //   else
  //   {
  //     sortKeyValuePair(_keys,_vals,_out_keys,_out_vals,len);

  //     for(unsigned int c_idx=C_row_ptr[i]; c_idx<C_row_ptr[i+1]; c_idx++)
  //     {
  //       C_col_ids[c_idx]  = _keys[c_idx - C_row_ptr[i]];
  //       C[c_idx]          = _vals[c_idx - C_row_ptr[i]];
  //     }
  //   }
  //   resetL1(_begg,_next,_keys,_vals,_pow,_size,_capacity);
  // }

}
