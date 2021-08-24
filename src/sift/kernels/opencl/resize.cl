__kernel void resize_nn(
    __global unsigned char *pixels,
    __global unsigned char *out,
    const int cols,
    const int rows,
    const int numOfChannels
)
{   
    int idx = get_global_id(0);
    int y = (idx/(numOfChannels)) / (cols);
    int x = (idx/(numOfChannels)) % (cols);
    int colorOffset = idx % numOfChannels;

    if (colorOffset != numOfChannels - 1) {
        out[idx] = (unsigned char)(y+x*(cols) + (y+2)+(x)*cols + y+(x+1)*(cols) + (y+2)+(x+1)*(cols))/4;
    } else {
        out[idx] = 255;
    }
}