#include <opencv2/core/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/features2d/features2d.hpp>

using namespace cv;

std::vector<DMatch> matchTwoImages(Mat d1, Mat d2) {
    Ptr<DescriptorMatcher> matcher = DescriptorMatcher::create(DescriptorMatcher::BRUTEFORCE);
    std::vector<DMatch> matches;
    matcher->match(d1, d2, matches);

    return matches;
}

void drawMatches(const char* fileName, Mat img1, Mat img2, std::vector<DMatch> matches, std::vector<cv::KeyPoint> k1, std::vector<cv::KeyPoint> k2) {
    Mat output;
    drawMatches(img1, k1, img2, k2, matches, output);
    imwrite(fileName, output);
}