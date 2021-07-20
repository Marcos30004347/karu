#pragma once

#include <string>

namespace karu {
namespace image {

class Image {
    public:
        Image(std::string filename, int num_of_channels);
        ~Image();

    private:
        int width;
        int height;
        int comp_per_pixel;
        unsigned char *data;
};


}
}