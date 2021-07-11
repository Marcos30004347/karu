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
	
  __global float *A,
	__global float *B,
	__global float *C,

  __local float* shared_x
) {
  const unsigned int t_block_idx = get_group_id(0);
  const unsigned int t_block_dim = get_local_size(0);
  const unsigned int t_idx      = get_local_id(0);

  const idx = t_block_idx * t_block_dim + t_idx;

  const unsigned long int bs =  block_heigth;

  unsigned long int A_target_block_row = t_block_idx/bs;
  unsigned long int B_target_block_column = t_block_idx/bs;
  unsigned long int r = t_idx;
  unsigned long int lane = t_idx;

  unsigned long int A_first_block = A_row_ptr[A_target_block_row];
  unsigned long int A_last_block = A_row_ptr[A_target_block_row+1];

  unsigned long int t = (idx / bs) % bs;

  // // block column
  // B[idx] = (idx/A_lines)%block_width;
  // // block row
  // B[idx] = t_block_idx / bs;

  // // block column
  // B[idx] = t_block_idx % bs;
  
  // block row
  // B[idx] = t_block_idx / bs;
  // B[idx] = lane/bs;
  // return;
  B[idx] = 0;
  unsigned long int block = A_first_block;
  // B[idx] = block*bs*bs;//B_col_ids[block];
  for(unsigned long int A_block = A_first_block; A_block < A_last_block; A_block++)
  {
    unsigned long int A_col = A_col_ids[A_block]/bs;
    unsigned long int A_lin = A_target_block_row;
  
    for(unsigned long int B_block = A_first_block; B_block < A_last_block; B_block++)
    {
      unsigned long int B_col = B_col_ids[A_block]/bs;
      
    }
    // Bs
    // Todo  
    unsigned long c = t_idx % bs;
    unsigned long l = t_idx / bs;
    B[A_block*bs*bs + c*bs + l] = (A_lin * 100 + A_col);
  }
}
