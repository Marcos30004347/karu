#include <bits/stdc++.h>
#include "Image.hpp"
#define STB_IMAGE_IMPLEMENTATION
#define STBI_FAILURE_USERMSG
#include "lib/stb_image.h"

using namespace karu;
using namespace image;

Image::Image(std::string filename, int num_of_channels) {
    int x,y,n;
    
    stbi_failure_reason();
    unsigned char *data = stbi_load(&filename[0], &x, &y, &n, num_of_channels);
    stbi_failure_reason();

    std::cout << x << " " << y << std::endl;

    assert(data != NULL);

    float *data_as_float;

    // convert char to float
    for (int i=0; i<y*x*num_of_channels; i++) {
        data_as_float[i] = (float)(data[i] - '0');
    }

    this->width = x;
    this->height = y;
    this->comp_per_pixel = n;
    this->data = data;
}

Image::~Image() {
    stbi_image_free(this->data);
}