#pragma once
#include "algebra/compute/Buffer.hpp"
#include "image/Image.hpp"

namespace karu {

class GaussianBlur {
private:
    float sigma;
    u64 numOfChannels;
    algebra::compute::Buffer *pixels;
    u64 rowsSize;
    u64 colSize;
    u64 kernelDimension;
    
public:
    GaussianBlur(float sigma, u64 numOfChannels, algebra::compute::Buffer *pixels, u64 rowSize, u64 colSize);
    GaussianBlur(float sigma, algebra::compute::Buffer *pixels, Image* img);
    ~GaussianBlur();
    void calculateKernel(f32 *GKernel);
    void run(algebra::compute::Buffer *out);

};

}
