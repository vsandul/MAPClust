#pragma once

#include <vector>
#include <set>
#include <map>
#include <string>

#include "../../hclust-cpp-master/fastcluster.h"
#include "ccl.h"
#include "dbscan.h"


struct SingleParticle{
    int Id;
    double X;
    double Y;
    double Z;
    double Px;
    double Py;
    double Pz;
    double E;
};

class BunchOfParticles{
public:
    BunchOfParticles(const std::vector<SingleParticle>& bunch):
        bunch_of_particles(bunch),
        labels_of_particles(bunch.size(), 0){}
    
    void ClusterizeBySpace(const double sp_dist, const std::string method );
    void ClusterizeByVelocity(const double v_dist, const std::string method );
    void Clusterize(const double sp_dist, const std::string sp_method, const double v_dist, const std::string v_method, const short clust_order = 0);
    
    const std::vector<int> GetLabelsVector();
    const std::map<int, std::vector<SingleParticle>> GetLabelsMap();    
    const std::vector<SingleParticle> GetBunchOfParticles();
    
private:
    const std::vector<SingleParticle> bunch_of_particles;
    std::vector<int> labels_of_particles;
    
    const double SpaceDist(const SingleParticle& part1, const SingleParticle& part2);
    const double VelDist(const SingleParticle& part1, const SingleParticle& part2);

    double** SpaceDistMat(const std::vector<SingleParticle>& bunch);
    double** VelDistMat(const std::vector<SingleParticle>& bunch);
    double* SpaceDistVec(const std::vector<SingleParticle>& bunch );
    double* VelDistVec(const std::vector<SingleParticle>& bunch);
    
    double** SpaceDistMat(const std::vector<int>& part_nums);
    double** VelDistMat(const std::vector<int>& part_nums);
    double* SpaceDistVec(const std::vector<int>& part_nums);
    double* VelDistVec(const std::vector<int>& part_nums);
    
    void ClusterizeBySpace(const double sp_dist, const std::string method, const std::vector<int>& part_orders, std::vector<int>& labels);
    void ClusterizeByVelocity(const double v_dist, const std::string method, const std::vector<int>& part_orders, std::vector<int>& labels);
    void ClusterizeBySpaceThanByVelocity(const double sp_dist, const std::string sp_method, const double v_dist, const std::string v_method);
    void ClusterizeByVelocityThanBySpace(const double sp_dist, const std::string sp_method, const double v_dist, const std::string v_method);
    void ClusterizeBySpaceAndVelocityTogether(const double sp_dist, const std::string sp_method, const double v_dist, const std::string v_method);
    
    const std::set<std::string> method_names = {
       "CCL", "DBSCAN", "HClust_Single", "HClust_Complete", "HClust_Average", "HClust_Median",        
    };
    
    const std::map<std::string, short> HClust_method_nummering = {
        {"HClust_Single", HCLUST_METHOD_SINGLE}, // = 0, single link with the minimum spanning tree algorithm (Rohlf, 1973)
        {"HClust_Complete", HCLUST_METHOD_COMPLETE}, // = 1, complete link with the nearest-neighbor-chain algorithm (Murtagh, 1984)
        {"HClust_Average", HCLUST_METHOD_AVERAGE}, // = 2, complete link with the nearest-neighbor-chain algorithm (Murtagh, 1984)
        {"HClust_Median", HCLUST_METHOD_MEDIAN}, // = 4, median link with the generic algorithm (Müllner, 2011)
    };

};
