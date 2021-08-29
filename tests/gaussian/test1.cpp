#include <string.h>
#include <iostream>

// #include "algebra/compute/Compute.hpp"
// #include "sift/GaussianBlur.hpp"
// #include "image/Image.hpp"

#include "sift/Sift.hpp"

using namespace karu;
// using namespace karu::algebra::compute;

int main() {

    //Context::initContext();

    Sift *sift = new Sift("../tests/assets/Karoo.png");

    sift->run();
    sift->drawImagesWithKeypoints("daleeee.jpg");

    return 0;
}