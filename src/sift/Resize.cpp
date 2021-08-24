#include <math.h>
#include "image/Image.hpp"
#include "sift/kernels/kernels.hpp"
#include "algebra/compute/Compute.hpp"

using namespace karu;
using namespace algebra::compute;

unsigned char* resizeImage(Image *img, f32 ratioX, f32 ratioY) {
    // Definir o tamanho da nova imagem
    int newWidth = ceilf(img->width/ratioX);
    int newHeight = ceilf(img->height/ratioY);
    int size = newWidth * newHeight * img->channels;

    Buffer pixels = img->createBuffer();
    Buffer out(size, Buffer::READ_WRITE, Buffer::MEM_GPU);

    pixels.upload();

    resize_kernel->setKernelArgument(0, pixels);
    resize_kernel->setKernelArgument(1, out);
    resize_kernel->setKernelArgument(2, sizeof(int), &img->width);
    resize_kernel->setKernelArgument(3, sizeof(int), &img->height);
    resize_kernel->setKernelArgument(4, sizeof(int), &img->channels);

    resize_kernel->enqueue({(u64)newHeight * newWidth * img->channels}, {1});
    
    unsigned char* newImg = (unsigned char*)out.download();
    return newImg;
}