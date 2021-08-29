#include <string.h>
#include <iostream>

// #include "algebra/compute/Compute.hpp"
// #include "sift/GaussianBlur.hpp"
// #include "image/Image.hpp"

#include "sift/Sift.hpp"
#include "sift/Matcher.cpp"

using namespace karu;
// using namespace karu::algebra::compute;

int main() {

    //Context::initContext();

    Sift *sift1 = new Sift("../tests/assets/Karoo.png");
    Sift *sift2 = new Sift("../tests/assets/Karoo.png");

    sift1->run();
    sift2->run();
    
    // sift1->drawImagesWithKeypoints("daleeee.jpg");
    // sift2->drawImagesWithKeypoints("daleeee.jpg");

    std::vector<DMatch> matches = matchTwoImages(sift1->src, sift2->src);
    drawMatches("daleeee.jpg",sift1->src, sift2->src, matches, sift1->keypoints, sift2->keypoints);

    return 0;
}