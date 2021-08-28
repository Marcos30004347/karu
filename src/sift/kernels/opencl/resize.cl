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

    printf("c: %u r: %u, channels: %u, idx: %i, y: %i, x: %i\n",cols,rows,numOfChannels,idx,y,x);

    if (colorOffset != numOfChannels - 1) {
        // out[idx] = (unsigned char)((pixels[y+x*(cols)] + pixels[(y+2)+(x)*cols] + pixels[y+(x+1)*(cols)] + pixels[(y+2)+(x+1)*(cols)])/4);
        out[idx] = pixels[idx];
    } else {
        out[idx] = 255;
    }
}
