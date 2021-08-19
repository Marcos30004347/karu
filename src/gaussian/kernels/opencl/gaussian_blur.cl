

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
    
    printf("colorOffset = %d\n", colorOffset);
    printf("cKernelDimension = %d\n", cKernelDimension);
    printf("currentRow = %d\n", currentRow);
    printf("currentCol = %d\n", currentCol);


    if (colorOffset != numOfChannels-1) {
        int i, j;
        
        for(j=0;j <=cKernelDimension; j++) {
            
            int y =currentRow + (j-(cKernelDimension/2));
            if(y < 0 || y >= rows) y = currentRow; 

            for(i=0;i <= cKernelDimension; i++) {
                
                int x = currentCol +(i-(cKernelDimension/2));
                if(x < 0 || x > cols) x = currentCol;
                acc += (float) ((float)(pixels[((y*(cols)+x))*numOfChannels + colorOffset])* ckernel[(j*(cKernelDimension))+i]) ;
            }
        }
        printf("%d, ", acc);
        if(acc >= 255) acc = 255;
        out[idx] = (unsigned char)acc;
    }  
    else {
        out[idx] = 255;
    }
    //out[idx] = pixels[idx];
}