#include "sift/Sift.hpp"

using namespace karu;

Sift::Sift(const char* filePath) {
    this->src = cv::imread(filePath, cv::IMREAD_GRAYSCALE);
    this->keypoints = {};
    this->descriptors = 0;
}

Sift::~Sift(){}

void Sift::run() {
    cv::Ptr<cv::SIFT> sift = cv::SIFT::create();
    sift->detectAndCompute(src,cv::noArray(),this->keypoints, this->descriptors);

}

void Sift::drawImagesWithKeypoints(const char* fileName) {
    cv::Mat output;
    cv::drawKeypoints(this->src, this->keypoints, output);
    cv::imwrite(fileName, output);
}