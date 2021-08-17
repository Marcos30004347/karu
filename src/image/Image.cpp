#include "image/Image.hpp"

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "libs/stb/stb_image.h"
#include "libs/stb/stb_image_write.h"

using namespace karu;
using namespace algebra::compute;

Image::Image(const char* filePath) {
    int width,height,channels;
    unsigned char *img = stbi_load(filePath, &width, &height, &channels, 0);
    
    this->pixels = img;
    this->width = width;
    this->height = height;
    this->channels = channels;
    this->imageSize = width * height * channels;

}

Image::Image(unsigned char* pixels, u64 width, u64 height, u64 channels) {
    this->pixels = pixels;
    this->width = width;
    this->height = height;
    this->channels = channels;
    this->imageSize = width * height * channels;
}

Image::Image(unsigned char* pixels, Image* img) {
    this->pixels = pixels;
    this->width = img->width;
    this->height = img->height;
    this->channels = img->channels;
    this->imageSize = this->width * this->height * this->channels;
}

Image::~Image() {
    free(this->pixels);
}

Image* Image::getGrayImage() {
    
    size_t img_size = this->width * this->height * this->channels;
    u64 gray_channels = this->channels == 4 ? 2 : 1;

    size_t gray_image_size = this->width * this->height * gray_channels;
    unsigned char *gray_img = new unsigned char[gray_image_size];

    if(gray_img == NULL) std::cout << "Erro alocação gray_img" << std::endl;

    for(unsigned char *p = this->pixels, *pg = gray_img; p != this->pixels + img_size; p += this->channels, pg += gray_channels) {
        *pg = (uint8_t)((*p + *(p+1) + *(p+2))/3.0);
        
        if(this->channels == 4) {
            *(pg+1) = *(p+3);
        }
    }

    return new Image(gray_img, this->width, this->height, gray_channels);
}

unsigned char* Image::getImage() {
    return this->pixels;
}

Buffer Image::createBuffer() {
    Buffer pixels(this->pixels, this->imageSize, Buffer::READ_WRITE, false);
    return pixels;
}

Buffer Image::createEmptyBuffer() {
    Buffer out(this->imageSize, Buffer::READ_WRITE, Buffer::MEM_GPU);
    return out;
}

void Image::writeImage(const char* path, int quality) {
    stbi_write_jpg(path, this->width, this->height, this->channels, this->pixels, quality);
}
