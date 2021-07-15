 __global__ void scan(
     float *g_odata,
     float *g_idata,
     int n
) {  
    extern __shared__ float temp[];
    
    int id = threadIdx.x;
    int pout = 0, pin = 1;

    temp[pout*n + id] = (id > 0) ? g_idata[id-1] : 0;
    
    __syncthreads();
    
    for (int offset = 1; offset < n; offset *= 2)
    {
        pout = 1 - pout;
        // swap double buffer indices
        pin = 1 - pout;

        if (id >= offset)
        {
            temp[pout*n + id] = temp[pin*n + id - offset] + temp[pout*n+id];
        }
        else 
        {
            temp[pout*n+id] = temp[pin*n+id];
            __syncthreads();
        }
        g_odata[id] = temp[pout*n+id];
        // write output 
    } 
}



__global__ void prescan(
    float *g_odata,
    float *g_idata,
    int n
) {
    extern __shared__ float temp[];

    int thid = threadIdx.x;
    int offset = 1; 
    
    temp[2*thid] = g_idata[2*thid];
    
    for (int d = n>>1; d > 0; d >>= 1)
    // build sum in place up the tree 
    { 
        __syncthreads();
        if (thid < d)
        { 
            int ai = offset*(2*thid+1)-1;
            int bi = offset*(2*thid+2)-1; 
            temp[bi] += temp[ai];
        }
        offset *= 2;
    }

    if (thid == 0)
    {
        temp[n - 1] = 0;
    }

    for (int d = 1; d < n; d *= 2)
    {
        offset >>= 1;
        __syncthreads();
        if (thid < d) 
        { 
            int ai = offset*(2*thid+1)-1;
            int bi = offset*(2*thid+2)-1; 
            float t = temp[ai];
            temp[ai] = temp[bi];
            temp[bi] += t;
        }
    }  
    __syncthreads(); 
    g_odata[2*thid] = temp[2*thid];
    // write results to device memory
    g_odata[2*thid+1] = temp[2*thid+1]; 
}

__kernel void prescan(
    __global int *g_odata,
    __global int *g_idata,
    __local int *temp,
    int n
) {
    int thid = get_local_id(0);
    int bid = get_group_id(0);
    int thread_num = get_local_size(0);

    int offset = 1;

    // Make the "empty" spots zeros, so it won't affect the final result.
    if((bid * thread_num + thid) < n)
    {
        temp[thid]  = g_idata[bid * thread_num + thid];
    }else
    {
        temp[thid]  = 0;
    } 

    // build sum in place up the tree
    for (int d = thread_num>>1; d > 0; d >>= 1)
    {
        barrier(CLK_LOCAL_MEM_FENCE);
        if (thid < d)
        {
            int ai = offset*(2*thid+1)-1;
            int bi = offset*(2*thid+2)-1;
            temp[bi] += temp[ai];
        }
        offset *= 2;
    }

    // clear the last element
    if (thid == 0)
    {
        temp[thread_num - 1] = 0;
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
            float t = temp[ai];
            temp[ai]  = temp[ bi];
            temp[bi] += t;
        }
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    g_odata[bid * thread_num + thid] = temp[thid];
}
