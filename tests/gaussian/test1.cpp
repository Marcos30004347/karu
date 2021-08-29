#include <string.h>
#include <iostream>

// #include "algebra/compute/Compute.hpp"
// #include "sift/GaussianBlur.hpp"
// #include "image/Image.hpp"

#include <opencv2/core/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/features2d/features2d.hpp>

// using namespace karu;
// using namespace karu::algebra::compute;

int main() {

    //Context::initContext();

    cv::Mat src = cv::imread("../tests/assets/Karoo.png", cv::IMREAD_GRAYSCALE);

    cv::Ptr<cv::SIFT> sift = cv::SIFT::create();
    std::vector<cv::KeyPoint> keypoints;
    sift->detect(src, keypoints);

    std::cout << keypoints[0] << std::endl;

    // Image *img = new Image("../tests/assets/Karoo.png");
    // Image *grayImg = img->getGrayImage();

    // Buffer pixels = grayImg->createBuffer();
    // Buffer out = grayImg->createEmptyBuffer();

    // pixels.upload();

    // GaussianBlur *gaussian = new GaussianBlur(5.0, &pixels, grayImg);
    // gaussian->run(&out);
    
    // unsigned char* new_img = (unsigned char*)out.download();
    // Image *newImg = new Image(new_img, grayImg); 

    // newImg->writeImage("luizin.jpg", 100);

    // Context::stopContext();

    return 0;
}