#include <math.h>
#include "image/Image.hpp"
#include "sift/kernels/kernels.hpp"
#include "algebra/compute/Compute.hpp"

using namespace karu;
using namespace algebra::compute;

unsigned char* resizeImage(Image *img, f32 ratioX, f32 ratioY) {
    // Definir o tamanho da nova imagem
    int newWidth = (int)ceil(img->width/ratioX);
    int newHeight = (int)ceil(img->height/ratioY);
    int size = newWidth * newHeight * (int)(img->channels);

    Buffer pixels = img->createBuffer();
    Buffer out(size*sizeof(unsigned char), Buffer::READ_WRITE, Buffer::MEM_GPU);

    pixels.upload();
    int largura = (int)img->width;
    int altura = (int)img->height;
    int canais = (int)img->channels;

    resize_kernel->setKernelArgument(0, pixels);
    resize_kernel->setKernelArgument(1, out);
    resize_kernel->setKernelArgument(2, sizeof(int), &largura);
    resize_kernel->setKernelArgument(3, sizeof(int), &altura);
    resize_kernel->setKernelArgument(4, sizeof(int), &canais);

    printf("=> %i\n",(newHeight * newWidth * (int)(img->channels)));
    resize_kernel->enqueue({(u64)(newHeight * newWidth * canais)}, {1});
    unsigned char* newImg = (unsigned char*)out.download();
    return newImg;
}
