#include "sift/GaussianBlur.hpp"
#include "sift/kernels/kernels.hpp"
#include "algebra/compute/Buffer.hpp"
#include <math.h>

#define debug(x) std::cout << #x << " = " << x << std::endl;

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

void GaussianBlur::calculateKernel(f32 *GKernel)
{
    i64 k = this->kernelDimension;
    double sigma = this->sigma;
    double r, s = 2.0 * sigma * sigma;
 
    // sum is for normalization
    double sum = 0.0;
 
    for (int x = -k/2; x <= k/2; x++) {
        for (int y = -k/2; y <= k/2; y++) {
            r = sqrt(x * x + y * y);
            // printf("%f ~ %f",(M_PI * s), (exp(-(r * r) / s)));
            GKernel[(x + k/2)*k + y + k/2] = (exp(-(r * r) / s)) / (M_PI * s);
            sum += GKernel[(x + k/2)*k + y + k/2];
        }
    }
 
    // normalising the Kernel
    for (int i = 0; i < k*k; ++i)
        GKernel[i] /= sum;
}


void GaussianBlur::run(Buffer *out) {

    u64 kernelDimension = (u64)ceilf(6 * sigma);
    if (kernelDimension % 2 == 0) kernelDimension++;

    this->kernelDimension = kernelDimension;

    f32 kernel[kernelDimension * kernelDimension];
    calculateKernel(kernel);
    Buffer new_kernel(kernel, sizeof(f32) * (this->kernelDimension * this->kernelDimension), Buffer::READ_WRITE, false);
    
    new_kernel.upload();

    // std::cout << "kernel dimension: " << kernelDimension << std::endl;
    // for (int i=0; i<kernelDimension; i++) {
    //     for (int j=0; j<kernelDimension; j++) {
    //         //std::cout << kernel[i] << ", " << std::endl;
    //         u64 aux =  i;
    //         u64 aux_2 = j;
    //         printf("%f ", kernel[aux * kernelDimension + aux_2]);
    //     }
    //     printf("\n");
    // }

    blur_kernel->setKernelArgument(0, this->pixels);
    
    blur_kernel->setKernelArgument(1, out);
    
    blur_kernel->setKernelArgument(2, new_kernel);

    blur_kernel->setKernelArgument(3, sizeof(u64), &this->rowsSize);
    
    blur_kernel->setKernelArgument(4, sizeof(u64), &this->colSize);
    
    blur_kernel->setKernelArgument(5, sizeof(u64), &this->kernelDimension);
    
    blur_kernel->setKernelArgument(6, sizeof(u64), &this->numOfChannels);
	
	  blur_kernel->enqueue({(u64)(this->rowsSize * this->colSize)}, {1});
    
    out->download();

}
