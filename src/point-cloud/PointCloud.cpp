#pragma once
#include "point-cloud/PointCloud.hpp"
#include "sift/Matcher.cpp"

using namespace karu;

PointCloud::PointCloud(std::vector<const char*> filesList) {
    this->filesList = filesList;
    this->sifts = {};
}

PointCloud::~PointCloud() {}

void PointCloud::run() {
    for (int i=0; i<filesList.size(); i++) {
        Sift *sift = new Sift(filesList[0]);
        sift->run();
        sifts.push_back(sift);
    }

    for (int i=0; i<this->sifts.size(); i++) {
        for (int j=i+1; j<this->sifts.size(); j++) {
            std::vector<cv::DMatch> matches = matchTwoImages(this->sifts[i]->descriptors, this->sifts[j]->descriptors);
            
        }
    }

}