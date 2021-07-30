
// 1 2 3 4   5 6 7 8
// 0       | 1
// 1 2 3 4 | 1 2 3 4
__kernel void reduce(
    __global float * input,
    __global float * output,
    __local float * target,
    unsigned int AOffset
) {
    const size_t globalId = get_global_id(0);
    const size_t localId  = get_local_id(0);
    target[localId] = input[globalId+AOffset];

    barrier(CLK_LOCAL_MEM_FENCE);

    size_t blockSize = get_local_size(0);
    size_t halfBlockSize = blockSize / 2;

    while (halfBlockSize>0) {
        if (localId<halfBlockSize) {
            target[localId] += target[localId + halfBlockSize];
            if ((halfBlockSize*2)<blockSize) { // uneven block division
                if (localId==0) { // when localID==0
                    target[localId] += target[localId + (blockSize-1)];
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        blockSize = halfBlockSize;
        halfBlockSize = blockSize / 2;
    }

    if (localId==0) {
        output[get_group_id(0)] = target[0];
    }
}
