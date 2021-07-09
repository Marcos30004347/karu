
unsigned long round_up_to_power_of_two(unsigned long v) {
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v++;
  return v;
}

__kernel void bsMV_kernel(
	const unsigned long int block_heigth,
	const unsigned long int block_width,
	const unsigned long int x_block_heigth,
	const unsigned long int x_block_width,
	__global unsigned long int *col_ids,
	__global unsigned long int *row_ptr,
	__global float *A,
	__global float *x,
	__global float *y,
  __local float* shared_x
) {
  const unsigned int t_block_idx = get_group_id(0);
  const unsigned int t_block_dim = get_local_size(0);
  const unsigned int t_idx = get_local_id(0);

  const idx = t_block_idx * t_block_dim + t_idx;

  const unsigned long int bs =  block_heigth;

  unsigned long int target_block_row = t_block_idx;
  unsigned long int r = t_idx;
  unsigned long int lane = t_idx;

  unsigned long int first_block = row_ptr[target_block_row];
  unsigned long int last_block = row_ptr[target_block_row+1];

  if(r < bs)
  {
    float local_out = 0;

    for(unsigned long int block = first_block; block < last_block; block++)
    {
      barrier(CLK_LOCAL_MEM_FENCE);
      shared_x[t_idx] = x[col_ids[block] + t_idx];
      barrier(CLK_LOCAL_MEM_FENCE);

      for(unsigned int c = 0; c<bs; c++)
      {
        local_out += shared_x[c] * A[block * bs*bs + c*bs + r];
      }
    }
  
    y[target_block_row * bs + r] = local_out;
  }
}
