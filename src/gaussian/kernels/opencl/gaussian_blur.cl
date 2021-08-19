

__kernel void blur(
    __global unsigned char *pixels,
    __global unsigned char *out,
    __global double* gaussian,
    const int rows,
    const int cols,
    const int k,
    const int numOfChannels
)
{
    int idx = get_global_id(0);
    int y = (idx/(numOfChannels))/ (cols);
    int x = (idx/(numOfChannels)) % (cols);
    int colorOffset = idx%(numOfChannels);
    float acc=0;
    
    // printf("colorOffset = %d\n", colorOffset);
    // printf("k = %d\n", k);
    // printf("y = %d\n", y);
    // printf("x = %d\n", x);
    // printf("num: %d\n",numOfChannels);

    
    int i, j;

    for (int c = 0; c < numOfChannels-1; c++) {
        // for(j=0;j <=k; j++) {
        
        //     int y =y + (j-(k/2));
        //     if(y < 0 || y >= rows) y = y; 

        //     for(i=0;i <= k; i++) {
                
        //         int x = x +(i-(k/2));
        //         if(x < 0 || x > cols) x = x;
        //         acc += (float) ((float)(pixels[((y*(cols)+x))*numOfChannels + colorOffset])* gaussian[(j*(k))+i]) ;
        //         //printf("acc: %d, ", acc);
        //     }
            
        // }

        int rowStart = max(0, y - k/2);
        int rowEnd = min(rows, y + k/2);

        int colStart = max(0, x - k/2);
        int colEnd = min(cols, x + k/2);

        int kernelColStart = min(x, k/2);
        int kernelColEnd = min(cols - x, k/2);

        int kernelRowStart = min(y, k/2);
        int kernelRowEnd = min(rows - y, k/2);

        //printf("%d, %d, %d, %d\n ", kernelRowStart, kernelRowEnd, kernelColStart, kernelColEnd);
        for (int row= k/2 - kernelRowStart; row < k/2 + kernelRowEnd; row++) {
            for (int col= k/2 - kernelColStart; col < k/2 + kernelColEnd; col++) {
                acc += gaussian[row*k+col] * pixels[((y-kernelRowStart + row) * cols + (x-kernelColStart + col)) * numOfChannels + c];
                //printf("%d, %d, %d\n ", row, col, gaussian[row*k+col]);
            }
        }
        //acc /= (k) * (k);
        out[(y * cols + x) * numOfChannels + c] = (unsigned char)acc;
    }

    out[(y * cols + x) * numOfChannels + numOfChannels - 1] = 255;

    //printf("%d, %d\n", pixels[(y * cols + x) * numOfChannels], pixels[(y * cols + x) * numOfChannels + 1]);
    //if(acc >= 255) acc = 255;
    //out[(y * cols + x) * numOfChannels + numOfChannels] = (unsigned char)255;

    //out[idx] = pixels[idx];
}