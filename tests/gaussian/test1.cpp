#include <string.h>
#include <iostream>

#include "algebra/compute/Compute.hpp"
#include "sift/GaussianBlur.hpp"
#include "image/Image.hpp"
#include "sift/Resize.cpp"

using namespace karu;
using namespace karu::algebra::compute;

int main() {

    Context::initContext();

    Image *img = new Image("../tests/assets/Karoo.png");
    Image *grayImg = img->getGrayImage();

    // unsigned char *pixelsss = grayImg->getImage();
    // for (int i=0; i<grayImg->imageSize; i++) {
    //     std::cout << (int)pixelsss[i] << ", ";
    // }


    // Buffer pixels = grayImg->createBuffer();
    // Buffer out = grayImg->createEmptyBuffer();

    // pixels.upload();

    // GaussianBlur *gaussian = new GaussianBlur(5.0, &pixels, grayImg);
    // gaussian->run(&out);
    // u64 newWidth;
    // u64 newHeight;

    u64 newWidth = ceilf(grayImg->width/2);
    u64 newHeight = ceilf(grayImg->height/2);

    unsigned char* new_img = resizeImage(grayImg, 2, 2);
    Image *newImg = new Image(new_img, newWidth, newHeight, grayImg->channels); 

    newImg->writeImage("luizin.jpg", 100);

    Context::stopContext();

    return 0;
}
