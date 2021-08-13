#include "gaussian/GaussianBlur.hpp"
#include "gaussian/kernels/kernels.hpp"
#include <math.h>

using namespace karu;

GaussianBlur::GaussianBlur(float sigma, int numOfChannels, unsigned char *pixels, int rowSize, int colSize) {
    this->sigma = sigma;
    this->numOfChannels = numOfChannels;
    this->pixels = pixels;
    this->rowsSize = rowSize;
    this->colSize = colSize;
}

GaussianBlur::~GaussianBlur(){}

float* GaussianBlur::calculateKernel() {
    float sigma = this->sigma;
    int kernelDimension = (int)ceilf(6 * sigma);
    
    if (kernelDimension % 2 == 0) kernelDimension++;
    int cKernelSize = pow(kernelDimension, 2);
    float *ckernel = new float[kernelDimension*kernelDimension];
 
    float acc = 0;
    
    for (int j = 0; j <= kernelDimension; j++) {
        int y = j - (kernelDimension / 2);
        
        for (int i = 0; i <= kernelDimension; i++) {
            int x = i - (kernelDimension / 2);
            ckernel[j*kernelDimension+i] = (1 / (2 * 3.14159*pow(sigma, 2)))*exp(-((pow(x, 2) + pow(y, 2)) / (2 * pow(sigma, 2)))); 
            acc += ckernel[j*i+i];
        }
    }
    
    for (int j = 0; j <= kernelDimension; j++) {
        for (int i = 0; i <= kernelDimension; i++) {
            ckernel[j*kernelDimension + i] = ckernel[j*kernelDimension + i] / acc;
        }
    }
    
    this->kernelDimension = kernelDimension;
}

unsigned char* GaussianBlur::run() {
    unsigned char* out;

    float* kernel = calculateKernel();
    
    blur_kernel->setKernelArgument(0, this->pixels);
    blur_kernel->setKernelArgument(1, out);
    blur_kernel->setKernelArgument(2, kernel);
    blur_kernel->setKernelArgument(3, this->rowsSize);
    blur_kernel->setKernelArgument(4, this->colSize);
    blur_kernel->setKernelArgument(5, this->kernelDimension);

    blur_kernel->enqueue({sizeof(this->pixels)}, {1});

    return out;
}