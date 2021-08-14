

__kernel void blur(
    __global unsigned char *pixels,
    __global unsigned char *out,
    __global float* ckernel,
    __constant int *rows,
    __constant int *cols,
    __constant int *cKernelDimension,
)
{
    int idx = get_global_id(0);
    int currentRow = (idx/(*channels))/ (*cols);
    int currentCol = (idx/(*channels)) % (*cols);
    int colorOffset = idx%(*channels);
    float acc=0;
 
    if (colorOffset != 3) {
        int i, j;
        
        for(j=0;j&amp;lt;(*cKernelDimension);j++) {
            
            int y =currentRow + (j-(*cKernelDimension/2));
            if(y < 0 || y >= *rows) y = currentRow; 

            for(i=0;i&amp;amp;lt;(*cKernelDimension);i++) {
                
                int x = currentCol +(i-(*cKernelDimension/2));
                if(x < 0 || x > *cols) x = currentCol;
                acc += (float) ((float)(pixels[((y*(*cols)+x))*4 + colorOffset])* ckernel[(j*(*cKernelDimension))+i]);
            }
        }
 
        if(acc &gt;= 255) acc = 255;
        out[idx] = (unsigned char)acc;
    }  
    else {
        out[idx] = 255;
    }
}