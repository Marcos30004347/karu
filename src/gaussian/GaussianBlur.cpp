#include "gaussian/GaussianBlur.hpp"
#include "gaussian/kernels/kernels.hpp"
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

f32 erfCustom(f32 x) {
    const f32 a1 = 0.254829592;
    const f32 a2 = -0.284496736;
    const f32 a3 = 1.421413741;
    const f32 a4 = -1.453152027;
    const f32 a5 = 1.061405429;
    const f32 p = 0.3275911;

    // A&S formula 7.1.26
    const f32 t = 1.0 / (1.0 + p * abs(x));
    const f32 y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x);

    debug(t);
    debug(y);

    if (x > 0) {
        x = 1;
    } else if (x < 0) {
        x = -1;
    } else {
        x = 0;
    }

    return x * y;
}

f32 defIntGaussian(f32 x, f32 mu, f32 sigma) {
    const f32 SQRT2 = 1.414;
    return 0.5 * erfCustom((x - mu) / (SQRT2 * sigma));
}

void gaussianKernel(f32 *coeff, u64 kernelSize, f32 sigma) {
    kernelSize = kernelSize*kernelSize;
    f32 mu = 0;
    f32 step = 1;
    
    const f32 end = 0.5 * kernelSize;
    const f32 start = -end;

    f32 sum = 0;
    f32 x = start;
    f32 lastInt = defIntGaussian(x, mu, sigma);
    f32 acc = 0;

    int idx = 0;

    while (x < end) {
        x += step;
        f32 newInt = defIntGaussian(x, mu, sigma);
        f32 c = newInt - lastInt;
        
        coeff[idx] = c;
        
        sum += c;
        lastInt = newInt;
        idx++;
    }

    // normalize
    sum = 1/sum;
    for (int i=0; i<idx+1; i++) {
        coeff[i] *= sum;
    }
}


void GaussianBlur::calculateKernel(f32 *GKernel) {
    // float sigma = this->sigma;
    
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

    u64 kernelDimension = 7;
    if (kernelDimension % 2 == 0) kernelDimension++;

    this->kernelDimension = kernelDimension;

    f32 kernel[kernelDimension * kernelDimension];
    gaussianKernel(kernel, kernelDimension, this->sigma);
    Buffer new_kernel(kernel, (u64)(this->kernelDimension * this->kernelDimension), Buffer::READ_WRITE, false);
    
    new_kernel.upload();

    std::cout << "kernel dimension: " << kernelDimension << std::endl;
    for (int i=0; i<kernelDimension*kernelDimension; i++) {
        std::cout << kernel[i] << ", " << std::endl;
    }

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