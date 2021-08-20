

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

    int i, j;

    for (int c = 0; c < numOfChannels-1; c++) {

        int rowStart = max(0, y - k/2);
        int rowEnd = min(rows, y + k/2);

        int colStart = max(0, x - k/2);
        int colEnd = min(cols, x + k/2);

        int kernelColStart = min(x, k/2);
        int kernelColEnd = min(cols - x, k/2);

        int kernelRowStart = min(y, k/2);
        int kernelRowEnd = min(rows - y, k/2);

        for (int row= k/2 - kernelRowStart; row < k/2 + kernelRowEnd; row++) {
            for (int col= k/2 - kernelColStart; col < k/2 + kernelColEnd; col++) {
                acc += gaussian[row*k+col] * pixels[((y-kernelRowStart + row) * cols + (x-kernelColStart + col)) * numOfChannels + c];
            }
        }
        out[(y * cols + x) * numOfChannels + c] = (unsigned char)acc;
    }

    out[(y * cols + x) * numOfChannels + numOfChannels - 1] = 255;

}