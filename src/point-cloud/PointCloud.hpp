#include "sift/Sift.hpp"
#include "bundle/BundleAdjustment.hpp"

namespace karu {

class PointCloud {
private:
    std::vector<const char*> filesList;
    std::vector<Sift*> sifts;

public:
    PointCloud(std::vector<const char*> filesList);
    ~PointCloud();
    void run();

};

}