#include <string.h>
#include <iostream>

#define STB_IMAGE_IMPLEMENTATION
#include "../gaussian/libs/stb/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../gaussian/libs/stb/stb_image_write.h"

#include "algebra/compute/Compute.hpp"
#include "gaussian/GaussianBlur.hpp"

using namespace karu;
using namespace karu::algebra::compute;

int main() {

    Context::initContext();

    int width,height,channels;
    // const char *filename = "Karoo.png";

    unsigned char *img = stbi_load("../tests/assets/Karoo.png", &width, &height, &channels, 0);
    
    std::cout << width << " " << height << " " << channels << std::endl;
    // std::cout << img << std::endl;

    // stbi_write_jpg("karu2.jpg",width,height,channels,img,1);

    size_t img_size = width * height * channels;
    u64 gray_channels = channels == 4 ? 2 : 1;

    size_t gray_image_size = width * height * gray_channels;
    unsigned char *gray_img = new unsigned char[gray_image_size];

    if(gray_img == NULL)
         std::cout << "Erro alocação gray_img" << std::endl;

    for(unsigned char *p = img, *pg = gray_img; p != img+img_size; p += channels, pg += gray_channels) {
        *pg = (uint8_t)((*p + *(p+1) + *(p+2))/3.0);
        if(channels == 4)
            *(pg+1) = *(p+3);
    }

    // stbi_write_jpg("karu_gray.jpg",width,height,gray_channels,gray_img,100);

    Buffer pixels(gray_img, gray_image_size, Buffer::READ_WRITE, false);
    Buffer out(gray_image_size, Buffer::READ_WRITE, Buffer::MEM_GPU);

    pixels.upload();

    GaussianBlur *gaussian = new GaussianBlur(2.0, gray_channels, &pixels, (u64)height, (u64)width);
    gaussian->run(&out);
    
    unsigned char* new_img = (unsigned char*)out.download();

    stbi_write_jpg("luizin.jpg",width,height,gray_channels, new_img,100);

    stbi_image_free(img);
    free(gray_img);
    free(new_img);

    Context::stopContext();

    return 0;
}