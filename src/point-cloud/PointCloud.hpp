#pragma once
#include "sift/Sift.hpp"
#include "algebra/compute/Compute.hpp"
#include <map>

namespace karu {

struct PointAux
{
    int img_idx;
    // {
    //     k1 -> {k3-img4},{k5-img2},...
    //     k2 -> ...
    //     ...
    // }
    std::map<u64,std::vector<std::pair<u64,u64>>> keyMap;
};

class PointCloud {
private:
    std::vector<const char*> filesList;
    std::vector<Sift*> sifts;
    std::vector<PointAux*> allMatches; 

public:
    PointCloud(std::vector<const char*> filesList);
    ~PointCloud();
    void run();

};

}