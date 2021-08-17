#pragma once
#include "algebra/compute/Compute.hpp"

namespace karu {

class Image {
private:
    unsigned char* pixels;

public:
    u64 width;
    u64 height;
    u64 channels;
    u64 imageSize;
    
    Image(const char* filePath);
    Image(unsigned char* pixels, u64 width, u64 height, u64 channels);
    Image(unsigned char* pixels, Image* img);

    ~Image();
    unsigned char* getImage();
    Image* getGrayImage();
    algebra::compute::Buffer createBuffer();
    algebra::compute::Buffer createEmptyBuffer();
    void writeImage(const char* path, int quality);
};

}