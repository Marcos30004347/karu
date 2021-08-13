#pragma once

namespace karu {

class GaussianBlur {
private:
    float sigma;
    int numOfChannels;
    unsigned char *pixels;
    int rowsSize;
    int colSize;
    int kernelDimension;
    
public:
    GaussianBlur(float sigma, int numOfChannels, unsigned char *pixels, int rowSize, int colSize);
    ~GaussianBlur();
    float* calculateKernel();
    unsigned char* run();

};

}
