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
        PointAux* aux = new PointAux;
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

    // u64 idx = 0;
    // std::vector<std::pair<u64,u64>> lista;
    // std::map<u64,std::vector<std::pair<u64,u64>>>::iterator it;
    // for (it=this->allMatches[0]->keyMap.begin(); it!=this->allMatches[0]->keyMap.end();  it++) {
    //     lista.push_back({it->first, idx});
    //     idx++;
    // }
    
    // para cada imagem
    // for (int i=1; i<this->sifts.size(); i++) {

    //     // para cada keymap dentro de cada PointAux
    //     for (it=this->allMatches[i]->keyMap.begin(); it!=this->allMatches[i]->keyMap.end();  it++) {
    //         bool flag = false;

    //         // para cada match do keypoint
    //         for(int j=0; j<it->second.size(); j++) {
    //             flag = false;

    //             // se deu match com alguma imagem anterior, colocar indice do keypoint "anterior"
    //             if(it->second[j].second < i) {
    //                 flag = true;
    //                 lista.push_back();
    //                 break;   
    //             }
    //         }

    //         // senao, keypoint novo, add na lista
    //         if(!flag){
    //             lista.push_back({it->first, idx});
    //             idx++;
    //         }
    //     }
    // }

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

}


/*
    [-1 -1 -1...] []
    [-1 -1 -1...]
    1 -> 2
    2 -> 3
    1 -> 3
*/

