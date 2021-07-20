#pragma once

#include "image/Image.hpp"

namespace karu {
namespace sift {



class Sift {
    public:
        void computeKeypointsAndDescriptors(image::Image *image);
    
    private:
        image::Image generateBaseImage(float sigma, float assumed_blur = 0.5);
};


}
}