
// #define SHARED_MEMORY_BANKS 32
// #define LOG_MEM_BANKS 5
// #define CONFLICT_FREE_OFFSET(n) ((n) >> LOG_MEM_BANKS)

int four_way_prefix_sum_with_shuffle_internal(
    __global int* keys_i,
    __global int* vals_i,
	__local int *cnt, 
	__local int *offsets,
	unsigned int blkSize,
	int thid,
	int gid,
	int extracted
)
{
	const unsigned int t_block_idx = get_group_id(0);
	
	/* Compute masks */
	for (int b = 0; b < 4; ++b)
	{
		cnt[b * blkSize + thid] = (extracted == b);
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	if (thid == 0)
	{
		for (int b = 0; b < 4; ++b)
		{
			offsets[b] = cnt[b * blkSize + blkSize - 1];
			// vals_i[t_block_idx*blkSize + thid] = cnt[thid * blkSize + blkSize - 1];
		}
	}

	int offset = 1;

	/*
	Prefix Sum in the Mask
	*/
	/* Build 4 ways sum tree */
	int n = blkSize;

	for (int d = n>>1; d > 0; d >>= 1)
	{
		barrier(CLK_LOCAL_MEM_FENCE);

		if (thid < d)
		{
			int ai = offset*(2*thid+1) - 1;
			int bi = offset*(2*thid+2) - 1;

			for (int b = 0; b < 4; ++b)
			{
				cnt[blkSize*b + bi] += cnt[blkSize*b + ai];
			}

		}

		offset*=2;
	}

	n = blkSize;

	/* Down-Sweep 4 ways */
	if (thid == 0)
	{
		for (int b = 0; b < 4; ++b)
			cnt[b*blkSize + n - 1] = 0;
	}

	for(int d=1; d < n; d*=2)
	{
		offset >>= 1;

		barrier(CLK_LOCAL_MEM_FENCE);

		if(thid < d)
		{
			int ai = offset*(2*thid+1) - 1;
			int bi = offset*(2*thid+2) - 1;
			
			barrier(CLK_LOCAL_MEM_FENCE);

			for(int b=0; b<4; b++)
			{
				int t = cnt[b*blkSize + ai];
				cnt[b*blkSize + ai] = cnt[b*blkSize + bi];
				cnt[b*blkSize + bi] += t;
			}
		}
	}

	barrier(CLK_LOCAL_MEM_FENCE);


	/* Get BlkSums */
	barrier(CLK_LOCAL_MEM_FENCE);
	
	if (thid == 0) {
		for (int b = 0; b < 4; ++b) {
			offsets[b] += cnt[b * blkSize + blkSize - 1];
		}
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	offset = 0;

	for (int i = 0; i < extracted; ++i)
		offset += offsets[i];

	return cnt[extracted * blkSize + thid] + offset;
}

__kernel void int_to_int_four_way_prefix_sum_shuffle(
	__global int* keys_i,
    __global int* vals_i,
    int bitIndx,
    int numElems,
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

	__local int offsets[4];

	// Do 4 way predication
	int key;
	int val;

	if (gid >= numElems && bitIndx == 0)
	{
		key = 2147483647;
		val = 0;
	}

	else
	{
		key = keys_i[gid];
		val = vals_i[gid];
	}

	int extracted = (key >> bitIndx) & 3;

	// Perform local shuffle and get local block sum data
	int addr = four_way_prefix_sum_with_shuffle_internal(keyShuffle_o, valShuffle_o, cnt, offsets, blkSize, thid, gid, extracted);

	s_keys[addr] = key;
	s_vals[addr] = val;

	keyShuffle_o[gid] = s_keys[thid];
	valShuffle_o[gid] = s_vals[thid];

	if (thid < 4) {
		blkSum_o[thid * get_num_groups(0) + get_group_id(0)] = offsets[thid];
	}

}


// int four_way_move_elements_internal
// (
// 	int thid, 
// 	__local int *counts, 
// 	__global int *blkSum, 
// 	__local int *prefixSums, 
// 	__global int *prefixBlkSum,
// 	__local int *offsets, 
// 	int extracted
// )
// {
// 	if (thid == 0) {
// 		offsets[0] = 0;

// 		for (int b = 0; b < 4; ++b) {
// 			counts[b] = blkSum[get_num_groups(0)*b + get_group_id(0)];
// 			prefixSums[b] = prefixBlkSum[get_num_groups(0)*b + get_group_id(0)];
// 			offsets[b+1] = offsets[b] + counts[b];
// 		}
// 	}

// 	barrier(CLK_LOCAL_MEM_FENCE);

// 	int Pdn = prefixSums[extracted];
// 	int m = thid - offsets[extracted];
// 	int a = Pdn + m;

// 	barrier(CLK_LOCAL_MEM_FENCE);

// 	return a;
// }

__kernel void move_int_to_int_elements(
	__global int *keyShuffle,
	__global int *valShuffle,
	__local int  *s_keys,
	__local int  *s_vals,
	__global int *blkSum,
	__global int *prefixBlkSum,
	__global int *keys,
	__global int *vals,
	int bitIndx,
	int numElems
)
{
	// Initialize
	int gid = get_global_id(0);
	int thid = get_local_id(0);

	__local int offsets[4];
	__local int counts[4];
	__local int prefixSums[4];

	// Get four way predication
	int key = keyShuffle[gid];
	int val = valShuffle[gid];

	int extracted = (key >> bitIndx) & 3;

	// Calculate the result address

	if (thid == 0) {
		offsets[0] = 0;

		for (int b = 0; b < 4; ++b) {
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

	if (pos < numElems) {
		keys[pos] = key;
		vals[pos] = val;
	}
}

// __kernel void parallel_order_checking(
// 	__global int *input,
// 	__global int *out,
// 	__local int  *temp,
// )
// {
// 	const unsigned int t_block_idx  = get_group_id(0);
// 	const unsigned int t_block_dim  = get_local_size(0);
// 	const unsigned int thid 		= get_local_id(0);

// 	unsigned int id = thid + t_block_idx*t_block_dim;

// 	temp[thid] = input[id];
// 	temp[thid+1] = input[id+1];

// 	barrier(CLK_LOCAL_MEM_FENCE);

// 	temp[thid] = temp[thid] > temp[thid+1];
// 	output[id] = temp[thid];
// }


// #define SHARED_MEMORY_BANKS 32
// #define LOG_MEM_BANKS 5

// // There were two BCAO optimisations in the paper - this one is fastest
// #define CONFLICT_FREE_OFFSET(n) ((n) >> LOG_MEM_BANKS)

// __kernel void stream_scan_kernel(
// 	__global int *output,
// 	__global int *input,
// 	__local int* temp,
// 	int n,
// 	int powerOfTwo
// ) {

// 	int threadID = get_global_id(0);
// 	// output[threadID] = threadID;

// 	// return;

// 	int ai = threadID;
// 	int bi = threadID + (n / 2);
// 	int bankOffsetA = CONFLICT_FREE_OFFSET(ai);
// 	int bankOffsetB = CONFLICT_FREE_OFFSET(bi);


// 	if (threadID < n) {
// 		temp[ai + bankOffsetA] = input[ai];
// 		temp[bi + bankOffsetB] = input[bi];
// 	}
// 	else {
// 		temp[ai + bankOffsetA] = 0;
// 		temp[bi + bankOffsetB] = 0;
// 	}


// 	int offset = 1;
// 	for (int d = powerOfTwo >> 1; d > 0; d >>= 1) // build sum in place up the tree
// 	{
// 		barrier(CLK_LOCAL_MEM_FENCE);
// 		if (threadID < d)
// 		{
// 			int ai = offset * (2 * threadID + 1) - 1;
// 			int bi = offset * (2 * threadID + 2) - 1;
// 			ai += CONFLICT_FREE_OFFSET(ai);
// 			bi += CONFLICT_FREE_OFFSET(bi);

// 			temp[bi] += temp[ai];
// 		}
// 		offset *= 2;
// 	}

// 	if (threadID == 0) {
// 		temp[powerOfTwo - 1 + CONFLICT_FREE_OFFSET(powerOfTwo - 1)] = 0; // clear the last element
// 	}

// 	for (int d = 1; d < powerOfTwo; d *= 2) // traverse down tree & build scan
// 	{
// 		offset >>= 1;
// 		barrier(CLK_LOCAL_MEM_FENCE);
// 		if (threadID < d)
// 		{
// 			int ai = offset * (2 * threadID + 1) - 1;
// 			int bi = offset * (2 * threadID + 2) - 1;
// 			ai += CONFLICT_FREE_OFFSET(ai);
// 			bi += CONFLICT_FREE_OFFSET(bi);

// 			int t = temp[ai];
// 			temp[ai] = temp[bi];
// 			temp[bi] += t;
// 		}
// 	}

// 	barrier(CLK_LOCAL_MEM_FENCE);

// 	if (threadID < n) {
// 		output[ai] = temp[ai + bankOffsetA];
// 		output[bi] = temp[bi + bankOffsetB];
// 	}
// }




// __kernel void _stream_scan_kernel(
//     __global int *g_odata,
//     __global int *g_idata,
// 		__global volatile int* intermediate,
//     __local int *temp,
//     int n
// ) {
// 	int sum = 0;
// 	int thid = get_local_id(0);
// 	int gid = get_global_id(0);
// 	int bid = get_group_id(0);
// 	int thread_num = get_local_size(0);

// 	g_odata[get_global_id(0)] = 0;
// 	// return;

// 	int offset = 1;

// 	// Make the empty spots zeros, so it won't affect the final result.
// 	if((bid * thread_num + thid) < n)
// 	{
// 			temp[thid]  = g_idata[bid * thread_num + thid];
// 	}else
// 	{
// 			temp[thid]  = 0;
// 	} 

// 	// build sum in place up the tree
// 	for (int d = thread_num>>1; d > 0; d >>= 1)
// 	{
// 			barrier(CLK_LOCAL_MEM_FENCE);
// 			if (thid < d)
// 			{
// 					int ai = offset*(2*thid+1)-1;
// 					int bi = offset*(2*thid+2)-1;
// 					temp[bi] += temp[ai];
// 			}
// 			offset *= 2;
// 	}

// 	barrier(CLK_LOCAL_MEM_FENCE);

// 	// clear the last element
// 	if(thid == 0)
// 	{
// 		intermediate[bid] = temp[thread_num - 1];
// 		temp[thread_num - 1] = 0;
// 	}

// 	// traverse down tree & build scan
// 	for (int d = 1; d < thread_num; d *= 2)
// 	{
// 			offset >>= 1;
// 			barrier(CLK_LOCAL_MEM_FENCE);
// 			if (thid < d)
// 			{
// 					int ai = offset*(2*thid+1)-1;
// 					int bi = offset*(2*thid+2)-1;
// 					float t = temp[ai];
// 					temp[ai]  = temp[ bi];
// 					temp[bi] += t;
// 			}
// 	}

// 	barrier(CLK_LOCAL_MEM_FENCE);

// 	g_odata[bid * thread_num + thid] = temp[thid];



// 	// if ((bid * thread_num + thid) < n) {
// 	// 	// Add sum of adjacent blocks
// 	// 	sum = (thid == thread_num - 1) ? temp[thread_num - 1] + temp[thread_num] : temp[thid];
	
// 	// 	if (bid != 0) sum += intermediate[bid - 1];
	
// 	// 	g_odata[(bid * thread_num + thid)] = sum;
// 	// }
	
// 	// const int B = (n / thread_num);

// 	// offset = 1;

// 	// for (int d = B>>1; d > 0; d >>= 1) // build sum in place up the tree
// 	// {
// 	// 		barrier(CLK_LOCAL_MEM_FENCE);
// 	// 		if (thid < d)
// 	// 		{
// 	// 			int ai = offset*(2*thid+1)-1;
// 	// 			int bi = offset*(2*thid+2)-1;
// 	// 			intermediate[bi] += intermediate[ai];
// 	// 		}
// 	// 		offset *= 2;
// 	// }

// 	// barrier(CLK_LOCAL_MEM_FENCE);

// 	// if (thid == 0 && bid == 0) {
// 	// 	intermediate[B - 1] = 0;
// 	// }

// 	// for (int d = 1; d < B; d *= 2) // traverse down tree & build scan
// 	// {
// 	// 	offset >>= 1;
// 	// 	barrier(CLK_LOCAL_MEM_FENCE);
// 	// 	if (thid < d)
// 	// 	{
// 	// 		int ai = offset*(2*thid+1)-1;
// 	// 		int bi = offset*(2*thid+2)-1;
// 	// 		int t = intermediate[ai];
// 	// 		intermediate[ai] = intermediate[bi];
// 	// 		intermediate[bi] += t;
// 	// 	}
// 	// }
// 	// 	barrier(CLK_LOCAL_MEM_FENCE);
	
// 	// // 	g_odata[2*thid] += aux[blockIdx.x];
// 	// // 	g_odata[2*thid+1] += aux[blockIdx.x];



// 	// g_odata[bid * thread_num + thid] = intermediate[bid];
// }



