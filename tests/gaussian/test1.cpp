#include <string.h>
#include <iostream>

#define STB_IMAGE_IMPLEMENTATION
#include "../gaussian/libs/stb/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../gaussian/libs/stb/stb_image_write.h"

#include "gaussian/GaussianBlur.hpp"

using namespace karu;

int main() {

    int width,height,channels;
    // const char *filename = "Karoo.png";
    unsigned char *img = stbi_load("./imgs/Karoo.png", &width, &height, &channels, 0);
    
    std::cout << width << " " << height << " " << channels << std::endl;
    // std::cout << img << std::endl;

    // stbi_write_jpg("karu2.jpg",width,height,channels,img,1);

    size_t img_size = width * height * channels;
    int gray_channels = channels == 4 ? 2 : 1;

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

    GaussianBlur *gaussian = new GaussianBlur(1.0, gray_channels, gray_img, height, width);
    unsigned char *new_img = gaussian->run();
    stbi_write_jpg("luizin.jpg",width,height,gray_channels,new_img,100);

    stbi_image_free(img);
    free(gray_img);
    free(new_img);
    
    // https://github.com/dbarac/sift-cpp
    return 0;
}