#include "point-cloud/PointCloud.hpp"
#include "bundle/BundleAdjustment.hpp"
#include "sift/Matcher.cpp"

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

    std::vector<Camera*> cameras;
    bool found = false;

    int img1 = -1;
    int img2 = -1;
    std::vector<std::pair<u64, u64>> pares;

    for (int i=0; i<this->allMatches.size(); i++) {
        
        //achar o indice de uma imagem que possui
        //pelo menos 8 matches com a imagem atual
        std::map<u64,std::vector<std::pair<u64,u64>>> keyMap = allMatches[i]->keyMap;
        std::map<u64,std::vector<std::pair<u64,u64>>>::iterator it;
        
        for (int j=0; j<this->sifts.size(); j++) {
            if (i == j) continue;
            int counter = 0;
            std::vector<std::pair<u64, u64>> curr_pares;
            
            for (it = keyMap.begin(); it != keyMap.end(); it++) {
                u64 keypoint = it->first;
                std::vector<std::pair<u64,u64>> keypointMatches = it->second;

                
                for (int k=0; k<keypointMatches.size(); k++) {
                    u64 keypointOfImg = keypointMatches[k].first;
                    u64 imgIdx = keypointMatches[k].second;
                    curr_pares.push_back(std::make_pair(keypoint, keypointOfImg));
                    if (imgIdx == j) counter++;
                    if (counter == 8) break;
                }
            }

            if (counter >= 8) {
                found = true;

                img1 = i;
                img2 = j;
                
                pares = curr_pares;

                break;
            }
        }

        if (found) break;

    }

    std::vector<Matrix> projections1, projections2;
    
    for (int i=0; i<pares.size(); i++) {
        int key1 = pares[i].first;
        int key2 = pares[i].second;

        Sift *sift1 = this->sifts[img1];
        Sift *sift2 = this->sifts[img2];

        projections1.push_back(Matrix(3,1, {
            sift1->keypoints[key1].pt.x,
            sift1->keypoints[key1].pt.y,
            1
        }));

        projections2.push_back(Matrix(3,1, {
            sift2->keypoints[key2].pt.x,
            sift2->keypoints[key2].pt.y,
            1
        }));
    }
    
    std::vector<u64> pointIdx1 = globalKeypoints[img1];
    std::vector<u64> pointIdx2 = globalKeypoints[img2];

    Camera cam1 = Camera(
        0,
        0,
        this->sifts[img1]->src.size[0]/2,
        this->sifts[img1]->src.size[1]/2,
        Matrix(3,1, {0, 0, 0}),
        algebra::rotationMaxtrixToAxisAngle(getRotationMatrix(2*PI, 2*PI, 2*PI))
    );

    Camera cam2 = Camera(
        0,
        0,
        this->sifts[img2]->src.size[0]/2,
        this->sifts[img2]->src.size[1]/2,
        Matrix(3,1, {0, 0, 0}),
        algebra::rotationMaxtrixToAxisAngle(getRotationMatrix(2*PI, 2*PI, 2*PI))
    );
    
    std::vector<bundle::Bundle> bundles;

    // Add to bundles   
    bundles.push_back({
        cam1, 
        projections1,
        pointIdx1
    });

    bundles.push_back({
        cam2,
        projections2,
        pointIdx2
    });

    std::vector<std::vector<u64>> point_idx_to_camera(pares.size());
    std::vector<Matrix> wPoints;

	for(i64 j=0; j<bundles.size(); j++)
		for(u64 i : bundles[j].point_idx)
			point_idx_to_camera[i].push_back(j);

	bundleAdjustment(bundles, wPoints, point_idx_to_camera, 0.000001);

	for(u64 j=0; j<bundles.size(); j++)
	{
		for(u64 i=0; i<bundles[j].projections.size(); i++)
		{
			Matrix r = bundles[j].projections[i] - bundles[j].camera.projection(wPoints[bundles[j].point_idx[i]]);
			assert(norm(r) < 0.000001);
		}
	}

}