

__kernel void blur(
    __global unsigned char *pixels,
    __global unsigned char *out,
    __global float* ckernel,
    const int rows,
    const int cols,
    const int cKernelDimension,
    const int numOfChannels
)
{
    int idx = get_global_id(0);
    int currentRow = (idx/(numOfChannels))/ (cols);
    int currentCol = (idx/(numOfChannels)) % (cols);
    int colorOffset = idx%(numOfChannels);
    float acc=0;
    
    // printf("colorOffset = %d\n", colorOffset);
    // printf("cKernelDimension = %d\n", cKernelDimension);
    // printf("currentRow = %d\n", currentRow);
    // printf("currentCol = %d\n", currentCol);
    // printf("num: %d\n",numOfChannels);

    
    int i, j;

    for (int c = 0; c < numOfChannels-1; c++) {
        // for(j=0;j <=cKernelDimension; j++) {
        
        //     int y =currentRow + (j-(cKernelDimension/2));
        //     if(y < 0 || y >= rows) y = currentRow; 

        //     for(i=0;i <= cKernelDimension; i++) {
                
        //         int x = currentCol +(i-(cKernelDimension/2));
        //         if(x < 0 || x > cols) x = currentCol;
        //         acc += (float) ((float)(pixels[((y*(cols)+x))*numOfChannels + colorOffset])* ckernel[(j*(cKernelDimension))+i]) ;
        //         //printf("acc: %d, ", acc);
        //     }
            
        // }
        
        out[(currentRow * cols + currentCol) * numOfChannels + c] = pixels[(currentRow * cols + currentCol) * numOfChannels + c];
    }

    out[(currentRow * cols + currentCol) * numOfChannels + numOfChannels - 1] = 255;

    //printf("%d, %d\n", pixels[(currentRow * cols + currentCol) * numOfChannels], pixels[(currentRow * cols + currentCol) * numOfChannels + 1]);
    //if(acc >= 255) acc = 255;
    //out[(currentRow * cols + currentCol) * numOfChannels + numOfChannels] = (unsigned char)255;

    //out[idx] = pixels[idx];
}