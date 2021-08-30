#pragma once
#include "point-cloud/PointCloud.hpp"
#include "sift/Matcher.cpp"
#include <map>

using namespace karu;



PointCloud::PointCloud(std::vector<const char*> filesList) {
    this->filesList = filesList;
    this->sifts = {};
}

PointCloud::~PointCloud() {}

void PointCloud::run() {
    for (int i=0; i<filesList.size(); i++) {
        Sift *sift = new Sift(filesList[i]);
        sift->run();
        sifts.push_back(sift);
    }

    for (int i=0; i<this->sifts.size(); i++) {
        PointAux* aux = new PointAux();
        aux->img_idx = i;
        for (int j=i+1; j<this->sifts.size(); j++) {
            std::vector<cv::DMatch> matches = matchTwoImages(this->sifts[i]->descriptors, this->sifts[j]->descriptors);
            for(int k=0; k<matches.size(); k++) {
                std::pair<u64,u64> p = {matches[k].trainIdx, j};
                aux->keyMap[matches[k].queryIdx].push_back(p);
            }
        }
        this->allMatches.push_back(aux);
    }

    std::vector<std::vector<u64>> globalKeypoints;

    // Passando por cada Point Aux e criando um vetor com o tamanho de todos os keypoints que deram match
    for (int i=0; i<this->allMatches.size(); i++) {
        std::map<u64,std::vector<std::pair<u64,u64>>> keyMap = allMatches[i]->keyMap;
        std::vector<u64>  vect(keyMap.size(), -1);
        globalKeypoints.push_back(vect);
    }

    u64 idx = 0;
    for (int i=0; i<this->allMatches.size(); i++) {
        std::map<u64,std::vector<std::pair<u64,u64>>> keyMap = allMatches[i]->keyMap;
        std::map<u64,std::vector<std::pair<u64,u64>>>::iterator it;
        
        for (it = keyMap.begin(); it != keyMap.end(); it++) {
            u64 keypoint = it->first;
            std::vector<std::pair<u64,u64>> keypointMatches = it->second;

            if (globalKeypoints[i][keypoint] == -1) {
                globalKeypoints[i][keypoint] = idx;
                idx++;
            }

            for (int m=0; m<keypointMatches.size(); m++) {
                u64 keypointOfImg = keypointMatches[m].first;
                u64 ImgIdx = keypointMatches[m].second;

                if (globalKeypoints[ImgIdx][keypointOfImg] == -1) {
                    globalKeypoints[ImgIdx][keypointOfImg] = globalKeypoints[i][keypoint];
                }
            }
        }
    }

    std::vector<bundle::Bundle> bundles;

    // Create Bundles
    for (int i=0; i<this->sifts.size(); i++) {
        // Criar a camera


        // Criar as projections
        Sift *sift = this->sifts[i];
        std::vector<Matrix> projections;
        for (const cv::KeyPoint& keypoint : sift->keypoints) {
            float x = keypoint.pt.x;
            float y = keypoint.pt.y;
            projections.push_back(Matrix(3,1, { x, y, 1 }));
        }

        // Add o global keypoints da imagem
        std::vector<u64> point_idx = globalKeypoints[i];

        // Add to bundles   
        bundles.push_back({

        });
    }

}


/*
    [-1 -1 -1...] []
    [-1 -1 -1...]
    1 -> 2
    2 -> 3
    1 -> 3
*/

