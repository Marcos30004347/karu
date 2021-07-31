__kernel void bsMV_kernel(
	const unsigned long int block_heigth,
	const unsigned long int block_width,
	const unsigned int x_block_heigth,
	const unsigned int x_block_width,
	const unsigned int x_stored_lines,
	const unsigned int x_stored_column,
	const unsigned int y_block_heigth,
	const unsigned int y_block_width,
	const __global unsigned long int *col_ids,
	const __global unsigned long int *row_ptr,
	const __global float *A,
	const __global float *x,
	__global float *y,
  __local float* shared_x
) {
  const unsigned int t_block_idx = get_group_id(0);
  const unsigned int t_block_dim = get_local_size(0);
  const unsigned int t_idx       = get_local_id(0);

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
			unsigned long int pos = col_ids[block] + t_idx;
      
			unsigned long int block_y = pos/x_block_heigth;
			unsigned long int by      = pos - block_y;
  
			unsigned long int blk_line_cnt = x_stored_column/x_block_width;
			unsigned long int x_blk_size   = x_block_heigth*x_block_width;
			unsigned long int stride       = block_y*blk_line_cnt*x_blk_size;
      

      barrier(CLK_LOCAL_MEM_FENCE);
      shared_x[t_idx] = x[stride + by*x_block_width];
      barrier(CLK_LOCAL_MEM_FENCE);
      

      for(unsigned int c = 0; c<bs; c++)
      {
        local_out += shared_x[c] * A[block*bs*bs + c*bs + r];
      }
    }

    unsigned long int pos = target_block_row + t_idx;
    unsigned long int block_y = pos/x_block_heigth;
    unsigned long int by      = pos - block_y*x_block_heigth;
    unsigned long int block_x = 0;
    unsigned long int bx      = 0 - block_x*x_block_width;
    unsigned long int blocks_per_block_line = x_stored_column/x_block_width;
    unsigned long int x_block_size = x_block_heigth * x_block_width;
    unsigned long int stride = block_y*blocks_per_block_line*x_block_size + block_x*x_block_size;
    
    y[target_block_row * bs + r] = local_out;

  }
}
