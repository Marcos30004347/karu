#include <string.h>
#include <iostream>

#include "algebra/compute/Compute.hpp"
#include "sift/GaussianBlur.hpp"
#include "image/Image.hpp"

using namespace karu;
using namespace karu::algebra::compute;

int main() {

    Context::initContext();

    Image *img = new Image("../tests/assets/Karoo.png");
    Image *grayImg = img->getGrayImage();

    Buffer pixels = grayImg->createBuffer();
    Buffer out = grayImg->createEmptyBuffer();

    pixels.upload();

    GaussianBlur *gaussian = new GaussianBlur(5.0, &pixels, grayImg);
    gaussian->run(&out);
    
    unsigned char* new_img = (unsigned char*)out.download();
    Image *newImg = new Image(new_img, grayImg); 

    newImg->writeImage("luizin.jpg", 100);

    Context::stopContext();

    return 0;
}
