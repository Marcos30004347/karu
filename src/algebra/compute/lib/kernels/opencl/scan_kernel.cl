// #define SHARED_MEMORY_BANKS 32
// #define LOG_MEM_BANKS 5
// #define CONFLICT_FREE_OFFSET(n) ((n) >> LOG_MEM_BANKS)

// This function scan an entire array, it is assume that 
// the length of the array is small than the block size
__kernel void small_scan_kernel(
	__global int *input,
	__global int *output,
	__local int* temp,
	int n,
	int powerOfTwo
) {
	int thid = get_local_id(0);

	// int ai = thid;
	// int bi = thid + (n / 2);

	// int bankOffsetA = CONFLICT_FREE_OFFSET(ai);
	// int bankOffsetB = CONFLICT_FREE_OFFSET(bi);

	if(thid < n)
	{
		temp[thid]  = input[thid];
	}else
	{
		temp[thid]  = 0;
	}

	int offset = 1;
	for (int d = powerOfTwo >> 1; d > 0; d >>= 1)
	{
		barrier(CLK_LOCAL_MEM_FENCE);
		if (thid < d)
		{
			int ai = offset * (2 * thid + 1) - 1;
			int bi = offset * (2 * thid + 2) - 1;

			// ai += CONFLICT_FREE_OFFSET(ai);
			// bi += CONFLICT_FREE_OFFSET(bi);

			temp[bi] += temp[ai];
		}
		offset *= 2;
	}

	if (thid == 0) {
		temp[powerOfTwo - 1 /*+ CONFLICT_FREE_OFFSET(powerOfTwo - 1)*/] = 0; // clear the last element
	}

	for (int d = 1; d < powerOfTwo; d *= 2)
	{
		offset >>= 1;
		barrier(CLK_LOCAL_MEM_FENCE);
		if (thid < d)
		{
			int ai = offset * (2 * thid + 1) - 1;
			int bi = offset * (2 * thid + 2) - 1;
			// ai += CONFLICT_FREE_OFFSET(ai);
			// bi += CONFLICT_FREE_OFFSET(bi);
			int t = temp[ai];
			temp[ai] = temp[bi];
			temp[bi] += t;
		}
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	if (thid < n) {
		output[thid] = temp[thid];
	}
}


__kernel void block_scan_kernel(
    __global int *g_idata,
    __global int *g_odata,
    __global volatile int* sums,
    __local int *temp,
    int n
) {
	int thid = get_local_id(0);
	int gid = get_global_id(0);
	int bid = get_group_id(0);
	int thread_num = get_local_size(0);

	// The worksize may be grather than the total array size,
	// so we only load to temp the avaliable values, 
	// 0 is loaded otherwise
	if(gid < n)
	{
		temp[thid]  = g_idata[gid];
	}else
	{
		temp[thid]  = 0;
	} 

	int offset = 1;

	// build sum in place up the tree
	for (int d = thread_num>>1; d > 0; d >>= 1)
	{
		barrier(CLK_LOCAL_MEM_FENCE);
		if (thid < d)
		{
			int ai = offset*(2*thid+1)-1;
			int bi = offset*(2*thid+2)-1;

			// ai += CONFLICT_FREE_OFFSET(ai);
			// bi += CONFLICT_FREE_OFFSET(bi);

			temp[bi] += temp[ai];
		}
		offset *= 2;
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	// clear the last element
	if(thid == 0)
	{
		sums[bid] = temp[thread_num - 1];
		temp[thread_num - 1 /* + CONFLICT_FREE_OFFSET(thread_num - 1) */] = 0;
	}

	// traverse down tree & build scan
	for (int d = 1; d < thread_num; d *= 2)
	{
		offset >>= 1;
		barrier(CLK_LOCAL_MEM_FENCE);
		if (thid < d)
		{
			int ai = offset*(2*thid+1)-1;
			int bi = offset*(2*thid+2)-1;
			// ai += CONFLICT_FREE_OFFSET(ai);
			// bi += CONFLICT_FREE_OFFSET(bi);

			int t = temp[ai];
			temp[ai]  = temp[ bi];
			temp[bi] += t;
		}
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	if(gid < n)
	{
		g_odata[gid] = temp[thid];
		// temp[thid]  = g_idata[gid];
	}

	return;
}


__kernel void sum_values(
    __global int *g_idata,
    __global int *g_odata,
    __global volatile int* sums
) {
	int thid = get_local_id(0);
	int gid = get_global_id(0);
	int bid = get_group_id(0);
	
	int thread_num = get_local_size(0);

	g_odata[gid] = g_idata[gid] + sums[bid];
}
