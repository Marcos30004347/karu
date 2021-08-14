#include "gaussian/GaussianBlur.hpp"
#include "gaussian/kernels/kernels.hpp"
#include "algebra/compute/Buffer.hpp"
#include <math.h>

using namespace karu;
using namespace algebra::compute;

GaussianBlur::GaussianBlur(float sigma, u64 numOfChannels, algebra::compute::Buffer *pixels, u64 rowSize, u64 colSize) {
    this->sigma = sigma;
    this->numOfChannels = numOfChannels;
    this->pixels = pixels;
    this->rowsSize = rowSize;
    this->colSize = colSize;
}

GaussianBlur::GaussianBlur(float sigma, Buffer *pixels, Image* img) {
    this->sigma = sigma;
    this->numOfChannels = img->channels;
    this->pixels = pixels;
    this->rowsSize = img->height;
    this->colSize = img->width;
}

GaussianBlur::~GaussianBlur(){}

void GaussianBlur::calculateKernel(f32 *GKernel) {
    // float sigma = this->sigma;
    // u64 kernelDimension = (u64)ceilf(6 * sigma);
    
    // if (kernelDimension % 2 == 0) kernelDimension++;
    // int cKernelSize = pow(kernelDimension, 2);
    // float *ckernel = new float[kernelDimension*kernelDimension];
 
    // float acc = 0;
    
    // for (int j = 0; j <= kernelDimension; j++) {
    //     int y = j - (kernelDimension / 2);
        
    //     for (int i = 0; i <= kernelDimension; i++) {
    //         int x = i - (kernelDimension / 2);
    //         ckernel[j*kernelDimension+i] = (1 / (2 * 3.14159*pow(sigma, 2)))*exp(-((pow(x, 2) + pow(y, 2)) / (2 * pow(sigma, 2)))); 
    //         acc += ckernel[j*i+i];
    //     }
    // }
    
    // for (int j = 0; j <= kernelDimension; j++) {
    //     for (int i = 0; i <= kernelDimension; i++) {
    //         ckernel[j*kernelDimension + i] = ckernel[j*kernelDimension + i] / acc;
    //     }
    // }
    

    // return ckernel;

    
    this->kernelDimension = 5;
    // initialising standard deviation to 1.0
    double sigma = this->sigma;
    double r, s = 2.0 * sigma * sigma;
 
    // sum is for normalization
    double sum = 0.0;
 
    // generating 5x5 kernel
    for (int x = -2; x <= 2; x++) {
        for (int y = -2; y <= 2; y++) {
            r = sqrt(x * x + y * y);
            GKernel[(x + 2)*this->kernelDimension + y + 2] = (exp(-(r * r) / s)) / (M_PI * s);
            sum += GKernel[(x + 2)*this->kernelDimension + y + 2];
        }
    }
 
    // normalising the Kernel
    for (int i = 0; i < this->kernelDimension; ++i)
        for (int j = 0; j < this->kernelDimension; ++j)
            GKernel[i*this->kernelDimension + j] /= sum;

}

void GaussianBlur::run(Buffer *out) {


    std::cout << "INIT" << std::endl;
    
    f32 kernel[25];
    calculateKernel(kernel);
    Buffer new_kernel(kernel, (u64)(this->kernelDimension * this->kernelDimension), Buffer::READ_WRITE, false);
    new_kernel.upload();

    std::cout << "CALC" << std::endl;

    blur_kernel->setKernelArgument(0, this->pixels);
    std::cout << "A" << std::endl;
    
    blur_kernel->setKernelArgument(1, out);
    std::cout << "B" << std::endl;
    
    blur_kernel->setKernelArgument(2, new_kernel);
    std::cout << "C" << std::endl;

    blur_kernel->setKernelArgument(3, sizeof(int), &this->rowsSize);
    std::cout << "D" << std::endl;
    
    blur_kernel->setKernelArgument(4, sizeof(int), &this->colSize);
    std::cout << "E" << std::endl;
    
    blur_kernel->setKernelArgument(5, sizeof(int), &this->kernelDimension);
    std::cout << "F" << std::endl;
    
    blur_kernel->setKernelArgument(6, sizeof(int), &this->numOfChannels);
    std::cout << "G" << std::endl;

    blur_kernel->enqueue({(u64)this->rowsSize * this->colSize * this->numOfChannels}, {1});
    
    out->download();
    
    std::cout << "H" << std::endl;

}