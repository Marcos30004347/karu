int roundPowerOfTwo(int v)
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

__kernel void int_to_int_four_way_prefix_sum_shuffle(
	__global int* keys_i,
  __global int* vals_i,
  int bitIndx,
  int size,
	__local int *cnt,
  __local int *s_keys,
  __local int *s_vals,
	__global int *blkSum_o,
  __global int *keyShuffle_o,
  __global int *valShuffle_o
) {
	// Initialize
	int gid = get_global_id(0);
	int thid = get_local_id(0);
	int blkSize = get_local_size(0);
	int offset = 1;

	__local int offsets[4];

	// Do 4 way predication
	int key = 2147483647;
	int val = 2147483647;

	if (gid < size)
	{
		key = keys_i[gid];
		val = vals_i[gid];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	/* Compute masks */
	int n = roundPowerOfTwo(blkSize);

	int extracted = (key >> bitIndx) & 3;

	for(int b=0; b<4; ++b)
	{
		cnt[b * n + thid] = (extracted == b);
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	if (thid == 0)
	{
		for(int b=0; b<4; ++b)
		{
			offsets[b] = cnt[b * n + n - 1];
		}
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	// Up-Seep 4 ways
	for (int d = n>>1; d > 0; d >>= 1)
	{

		// barrier(CLK_LOCAL_MEM_FENCE);
		if (thid < d)
		{
			int ai = offset*(2*thid+1) - 1;
			int bi = offset*(2*thid+2) - 1;
			for(int b=0; b<4; ++b)
			{
				cnt[n*b + bi] += cnt[n*b + ai];
			}
		}
		barrier(CLK_LOCAL_MEM_FENCE);

		offset*=2;
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	// Down-Sweep 4 ways
	if (thid == 0)
	{
		for(int b=0; b<4; ++b)
		{
			cnt[b*n + n - 1] = 0;
		}
	}

	barrier(CLK_LOCAL_MEM_FENCE);
	
	for(int d=1; d<n; d*=2)
	{

		offset >>= 1;

		if(thid < d)
		{

			int ai = offset*(2*thid+1) - 1;
			int bi = offset*(2*thid+2) - 1;

			for(int b=0; b<4; ++b)
			{
				int t = cnt[b*n + ai];
				cnt[b*n + ai] = cnt[b*n + bi];
				cnt[b*n + bi] += t;
			}
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	// build Block Sums

	barrier(CLK_LOCAL_MEM_FENCE);

	if (thid == 0) {
		for(int b=0; b<4; ++b)
		{
			offsets[b] += cnt[b*n + n-1];
		}
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	offset = 0;

	for (int i = 0; i < extracted; ++i)
		offset += offsets[i];

	barrier(CLK_LOCAL_MEM_FENCE);
	
	int addr = cnt[extracted * n + thid] + offset;

	s_keys[addr] = key;
	s_vals[addr] = val;

	barrier(CLK_LOCAL_MEM_FENCE);

	keyShuffle_o[gid] = s_keys[thid];
	valShuffle_o[gid] = s_vals[thid];

	if (thid < 4)
	{
		blkSum_o[thid * get_num_groups(0) + get_group_id(0)] = offsets[thid];
	}
}

__kernel void move_int_to_int_elements(
	__global int *keyShuffle,
	__global int *valShuffle,
	__global int *blkSum,
	__global int *prefixBlkSum,
	__global int *keys,
	__global int *vals,
	int bitIndx,
	int size
)
{
	// Initialize
	int gid = get_global_id(0);
	int thid = get_local_id(0);

	__local int offsets[5];
	__local int counts[4];
	__local int prefixSums[4];

	int key = keyShuffle[gid];
	int val = valShuffle[gid];

	int extracted = (key >> bitIndx) & 3;

	if (thid == 0)
	{
		offsets[0] = 0;
	
		for (int b = 0; b < 4; ++b)
		{
			counts[b] = blkSum[get_num_groups(0)*b + get_group_id(0)];
			prefixSums[b] = prefixBlkSum[get_num_groups(0)*b + get_group_id(0)];
			offsets[b+1] = offsets[b] + counts[b];
		}
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	int Pdn = prefixSums[extracted];
	int m = thid - offsets[extracted];
	int pos = Pdn + m;

	barrier(CLK_LOCAL_MEM_FENCE);
	
	if (pos < size)
	{
		keys[pos] = key;
		vals[pos] = val;
	}
}

__kernel void parallel_block_order_checking(
	__global int *input,
	__global volatile int* output,
	__local int  *temp
)
{
	const unsigned int t_block_idx  = get_group_id(0);
	const unsigned int t_block_dim  = get_local_size(0);
	const unsigned int thid 		= get_local_id(0);

	unsigned int id = thid + t_block_idx*t_block_dim;

	temp[thid] = input[id];

	if(thid == 0)
	{
		if(t_block_idx < (get_num_groups(0) - 1))
		{
			temp[t_block_dim] = input[t_block_idx*t_block_dim + t_block_dim];
		}
		else
		{
			temp[t_block_dim] = 2147483647;
		}
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	temp[thid] = temp[thid] > temp[thid+1];

	int blockSize = get_local_size(0);
	int halfBlockSize = blockSize / 2;

	while (halfBlockSize > 0) {
		if (thid < halfBlockSize) {
			temp[thid] += temp[thid + halfBlockSize];

			if ((halfBlockSize*2)<blockSize) { // uneven block division
				if (thid==0) { // when localID==0
					temp[thid] += temp[thid + (blockSize-1)];
				}
			}

		}
		barrier(CLK_LOCAL_MEM_FENCE);

		blockSize = halfBlockSize;
		halfBlockSize = blockSize / 2;
	}

	if (thid==0) {
			output[get_group_id(0)] = temp[0];
	}
}
