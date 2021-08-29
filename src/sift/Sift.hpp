#pragma once
#include "algebra/compute/Compute.hpp"
#include <opencv2/core/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/features2d/features2d.hpp>

namespace karu {

class Sift {
private:
    cv::Mat src;
    std::vector<cv::KeyPoint> keypoints;
    cv::Mat descriptors;

public:
    
    Sift(const char* filePath);
    ~Sift();

    void run();
    void drawImagesWithKeypoints(const char* fileName); 
    
};

}