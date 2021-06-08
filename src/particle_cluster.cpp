#include <iterator>
#include <iostream>
#include <utility>

#include "particle_cluster.h"
#include "relativistic_kinematics.h"


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// Calculate space and velocity distances of two particles //////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const double BunchOfParticles::SpaceDist(const SingleParticle& part1, const SingleParticle& part2){
    const std::vector<double> r1 = {part1.X, part1.Y, part1.Z};
    const std::vector<double> r2 = {part2.X, part2.Y, part2.Z};
    const std::vector<double> w = {part1.Vx, part1.Vy, part1.Vz};
    const std::vector<double> r1_ = CoordinateInNewFrame(r1, w);
    const std::vector<double> r2_ = CoordinateInNewFrame(r2, w);
    const auto dX = (r1_.at(0)-r2_.at(0));
    const auto dY = (r1_.at(1)-r2_.at(1));
    const auto dZ = (r1_.at(2)-r2_.at(2));
    return sqrt(dX*dX+dY*dY+dZ*dZ);        
}

const double BunchOfParticles::VelDist(const SingleParticle& part1, const SingleParticle& part2){
    const std::vector<double> v1 = {part1.Vx, part1.Vy, part1.Vz};
    const std::vector<double> v2 = {part2.Vx, part2.Vy, part2.Vz};
    //const std::vector<double> w = {part1.Vx, part1.Vy, part1.Vz};
    //const std::vector<double> v1_ = VelocityInNewFrame(v1, v1);
    const std::vector<double> v1_ = {0,0,0};
    const std::vector<double> v2_ = VelocityInNewFrame(v2, v1);
    const auto dVx = (v1_.at(0)-v2_.at(0));
    const auto dVy = (v1_.at(1)-v2_.at(1));
    const auto dVz = (v1_.at(2)-v2_.at(2));
    return sqrt(dVx*dVx+dVy*dVy+dVz*dVz);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////// Calculate distances matrix/vector by vector of particles /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double** BunchOfParticles::SpaceDistMat(const std::vector<SingleParticle>& bunch){
    int npoints = bunch.size(); 
    double** distmat = new double* [npoints];
	for (int i = 0; i < npoints; i++)
		distmat[i] = new double [npoints];
    if (npoints == 0) {
        std::cout << "Empty Cluster => empty distances matrix." << std::endl;
        return distmat;
    }
    for (size_t i=0; i<npoints; i++) { 
        for (size_t j=i; j<npoints; j++) {
            if (j == i)
                distmat[i][i] = 0;
            else
                distmat[i][j] = distmat[j][i] = SpaceDist(bunch.at(i), bunch.at(j));
        }
    }
    return distmat;    
}

double** BunchOfParticles::VelDistMat(const std::vector<SingleParticle>& bunch){
    int npoints = bunch.size();
    double** distmat = new double* [npoints];
	for (int i = 0; i < npoints; i++)
		distmat[i] = new double [npoints];
    if (npoints == 0) {
        std::cout << "Empty Cluster => empty distances matrix." << std::endl;
        return distmat;
    }
    for (size_t i=0; i<npoints; i++) {
        for (size_t j=i; j<npoints; j++) {
            if (j == i)
                distmat[i][i] = 0;
            else
                distmat[i][j] = distmat[j][i] = VelDist(bunch.at(i), bunch.at(j));
        }
    }
    return distmat; 
}

double* BunchOfParticles::SpaceDistVec(const std::vector<SingleParticle>& bunch){
    size_t counter = 0;
    size_t npoints = bunch.size();
    double* distvec = new double[npoints*(npoints-1)/2];
    if (npoints == 0) {
        std::cout << "Empty Cluster => empty distances vector." << std::endl;
        return distvec;
    } 
    if (npoints == 1 ){
        std::cout << "1-particle Cluster => empty distances vector." << std::endl;
        return distvec;
    }
    for (size_t i=0; i<npoints; i++) {
        for (size_t j=i+1; j<npoints; j++) {
            distvec[counter] = SpaceDist(bunch.at(i), bunch.at(j));
            counter++;
        }
    }    
    return distvec;
}

double* BunchOfParticles::VelDistVec(const std::vector<SingleParticle>& bunch){
    size_t counter = 0;
    size_t npoints = bunch.size();
    double* distvec = new double[npoints*(npoints-1)/2];
    if (npoints == 0) {
        std::cout << "Empty Cluster => empty distances vector." << std::endl;
        return distvec;
    }
    if (npoints == 1) {
        std::cout << "1-particle Cluster => empty distances vector." << std::endl;
        return distvec;
    }    
    for (size_t i=0; i<npoints; i++) {
        for (size_t j=i+1; j<npoints; j++) {
            distvec[counter] = VelDist(bunch.at(i), bunch.at(j));
            counter++;
        }
    }
    return distvec;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////// Calculate distances matrix/vector by vector of orders of particles in 'bunch_of_particles' ///////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double** BunchOfParticles::SpaceDistMat(const std::vector<int>& part_orders){
    int npoints = part_orders.size();
    double** distmat = new double* [npoints];
	for (int i = 0; i < npoints; i++)
		distmat[i] = new double [npoints];
    if (npoints == 0) {
        std::cout << "Empty Cluster => empty distances matrix." << std::endl;
        return distmat;
    }
    for (size_t i=0; i<npoints; i++) {
        for (size_t j=i; j<npoints; j++) {
            if (j == i)
                distmat[i][i] = 0;
            else
                distmat[i][j] = distmat[j][i] = SpaceDist(bunch_of_particles.at(part_orders.at(i)),
                                                          bunch_of_particles.at(part_orders.at(j)) );
        }
    }
    return distmat;    
}

double** BunchOfParticles::VelDistMat(const std::vector<int>& part_orders){
    int npoints = part_orders.size();
    double** distmat = new double* [npoints];
	for (int i = 0; i < npoints; i++)
		distmat[i] = new double [npoints];
    if (npoints == 0) {
        std::cout << "Empty Cluster => empty distances matrix." << std::endl;
        return distmat;
    }
    for (size_t i=0; i<npoints; i++) {
        for (size_t j=i; j<npoints; j++) {
            if (j == i)
                distmat[i][i] = 0;
            else
                distmat[i][j] = distmat[j][i] = VelDist(bunch_of_particles.at(part_orders.at(i)),
                                                        bunch_of_particles.at(part_orders.at(j)) );
        }
    }
    return distmat; 
}

double* BunchOfParticles::SpaceDistVec(const std::vector<int>& part_orders){
    size_t counter = 0;
    size_t npoints = part_orders.size();
    double* distvec = new double[npoints*(npoints-1)/2];
    if (npoints == 0) {
        std::cout << "Empty Cluster => empty distances vector." << std::endl;
        return distvec;
    } 
    if (npoints == 1) {
        std::cout << "1-particle Cluster => empty distances vector." << std::endl;
        return distvec;
    }
    for (size_t i=0; i<npoints; i++) {
        for (size_t j=i+1; j<npoints; j++) {
            distvec[counter] = SpaceDist(bunch_of_particles.at(part_orders.at(i)),
                                            bunch_of_particles.at(part_orders.at(j)) );
            counter++;
        }
    }     
    return distvec;
}

double* BunchOfParticles::VelDistVec(const std::vector<int>& part_orders){
    size_t counter = 0;
    size_t npoints = part_orders.size();
    double* distvec = new double[npoints*(npoints-1)/2];
    if (npoints == 0) {
        std::cout << "Empty Cluster => empty distances vector." << std::endl;
        return distvec;
    } 
    if (npoints == 1) {
        std::cout << "1-particle Cluster => empty distances vector." << std::endl;
        return distvec;
    }
    for (size_t i=0; i<npoints; i++) {
        for (size_t j=i+1; j<npoints; j++) {
            distvec[counter] = VelDist(bunch_of_particles.at(part_orders.at(i)),
                                        bunch_of_particles.at(part_orders.at(j)) );
            counter++;
        }
    }
    return distvec;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Clusterize some particles by space distances///////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void BunchOfParticles::ClusterizeBySpace(const double sp_dist, const std::string method ){
    int nparts = bunch_of_particles.size();  
    if (nparts == 0) {
        std::cout << "No particles to clusterize! Exit with doing nothing."<< std::endl;
        return;
    }
    if (nparts == 1){
        labels_of_particles = {0};
        return;
    }
    
    if (method_names.find(method) == method_names.end()){
        std::cerr << "BunchOfParticles::ClusterizeBySpace: Error! Wrong method name.\nAllowed methods: ";
        int counter = 0;
        for (const auto& name:method_names){
            if (counter++<method_names.size()-1)
                std::cerr << name << ", ";
            else
                std::cerr << name << '\n';
        }            
        throw ;   
    }
    
    if (method == "CCL"){
        double** distmat = SpaceDistMat(bunch_of_particles);
        labels_of_particles = CCLClust(distmat, nparts, sp_dist);            
        return;
    } else if (method == "DBSCAN") {
        double** distmat = SpaceDistMat(bunch_of_particles);
        labels_of_particles = DBSCAN( distmat, nparts, sp_dist, 1);
        return;
    } else {
        auto opt_method = HClust_method_nummering.at(method);
        double* distvec = SpaceDistVec(bunch_of_particles);  
        int* merge = new int[2*(nparts-1)];
        double* height = new double[nparts-1];
        hclust_fast(nparts, distvec, opt_method, merge, height);                
        int* labels = &labels_of_particles[0];
        cutree_cdist(nparts, merge, height, sp_dist, labels);
        return;
    }                
}

void BunchOfParticles::ClusterizeBySpace(const double sp_dist, const std::string method, 
                                         const std::vector<int>& part_orders, std::vector<int>& labels){ 
    int nparts = part_orders.size();  
    if (nparts == 0) {
        std::cout << "No particles to clusterize! Exit with doing nothing."<< std::endl;
        return;
    }
    if (nparts == 1){
        labels = {0};
        return;
    }
    
    if (method_names.find(method) == method_names.end()){
        std::cerr << "BunchOfParticles::ClusterizeBySpace: Error! Wrong method name.\nAllowed methods: ";
        int counter = 0;
        for (const auto& name:method_names){
            if (counter++<method_names.size()-1)
                std::cerr << name << ", ";
            else
                std::cerr << name << '\n';
        }            
        throw ;   
    }  

    if (method == "CCL"){
        double** distmat = SpaceDistMat(part_orders);
        labels = CCLClust(distmat, nparts, sp_dist);
        return;
    } else if (method == "DBSCAN") {
        double** distmat = SpaceDistMat(part_orders);
        labels = DBSCAN( distmat, nparts, sp_dist, 1);
        return;
    } else {
        auto opt_method = HClust_method_nummering.at(method);
        double* distvec = SpaceDistVec(part_orders);  
        int* merge = new int[2*(nparts-1)];
        double* height = new double[nparts-1];
        hclust_fast(nparts, distvec, opt_method, merge, height);                
        int* labels_ = &labels[0];
        cutree_cdist(nparts, merge, height, sp_dist, labels_);
        return;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Clusterize some particles by velocity distances////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void BunchOfParticles::ClusterizeByVelocity(const double v_dist, const std::string method ){
    int nparts = bunch_of_particles.size();  
    if (nparts == 0) {
        std::cout << "No particles to clusterize! Exit with doing nothing."<< std::endl;
        return;
    }
    if (nparts == 1){
        labels_of_particles = {0};
        return;
    }
    
    if (method_names.find(method) == method_names.end()){
        std::cerr << "BunchOfParticles::ClusterizeByVelocity: Error! Wrong method name.\nAllowed methods: ";
        int counter = 0;
        for (const auto& name:method_names){
            if (counter++<method_names.size()-1)
                std::cerr << name << ", ";
            else
                std::cerr << name << '\n';
        }            
        throw ;   
    }
    
    if (method == "CCL"){
        double** distmat = VelDistMat(bunch_of_particles);
        labels_of_particles = CCLClust(distmat, nparts, v_dist);
        return;
    } else if (method == "DBSCAN") {
        double** distmat = VelDistMat(bunch_of_particles);
        labels_of_particles = DBSCAN( distmat, nparts, v_dist, 1);
        return;
    } else {
        auto opt_method = HClust_method_nummering.at(method);
        double* distvec = VelDistVec(bunch_of_particles);  
        int* merge = new int[2*(nparts-1)];
        double* height = new double[nparts-1];
        hclust_fast(nparts, distvec, opt_method, merge, height);                
        int* labels = &labels_of_particles[0];
        cutree_cdist(nparts, merge, height, v_dist, labels);
        return;
    }
}


void BunchOfParticles::ClusterizeByVelocity(const double v_dist, const std::string method, 
                                            const std::vector<int>& part_orders, std::vector<int>& labels){
    int nparts = part_orders.size();  
    if (nparts == 0) {
        std::cout << "No particles to clusterize! Exit with doing nothing."<< std::endl;
        return;
    }
    if (nparts == 1){
        labels = {0};
        return;
    }
    
    if (method_names.find(method) == method_names.end()){
        std::cerr << "BunchOfParticles::ClusterizeByVelocity: Error! Wrong method name.\nAllowed methods: ";
        int counter = 0;
        for (const auto& name:method_names){
            if (counter++<method_names.size()-1)
                std::cerr << name << ", ";
            else
                std::cerr << name << '\n';
        }            
        throw ;   
    }   

    if (method == "CCL"){
        double** distmat = VelDistMat(part_orders);
        labels = CCLClust(distmat, nparts, v_dist);
        return;
    } else if (method == "DBSCAN") {
        double** distmat = VelDistMat(part_orders);
        labels = DBSCAN( distmat, nparts, v_dist, 1);
        return;
    } else {
        auto opt_method = HClust_method_nummering.at(method);
        double* distvec = VelDistVec(part_orders); 
        int* merge = new int[2*(nparts-1)];
        double* height = new double[nparts-1];
        hclust_fast(nparts, distvec, opt_method, merge, height);                
        int* labels_ = &labels[0];
        cutree_cdist(nparts, merge, height, v_dist, labels_); 
        return;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// Clusterize by order/////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void BunchOfParticles::ClusterizeBySpaceThanByVelocity(const double sp_dist, const std::string sp_method,
                                                       const double v_dist, const std::string v_method){
    if (bunch_of_particles.size() == 0) {
        std::cout << "No particles to clusterize! Exit with doing nothing."<< std::endl;
        return;
    }
    if (bunch_of_particles.size() == 1){
        labels_of_particles = {0};
        return;
    }
    
    std::map< std::pair<int, int>, std::vector<int> > map_of_labels;
    
    std::map <int, std::vector<int>> map_of_space_labels;    
    ClusterizeBySpace(sp_dist, sp_method);    
    for (size_t i = 0; i < labels_of_particles.size(); i++)
        map_of_space_labels[labels_of_particles.at(i)].push_back(i);
    
    for( auto it = map_of_space_labels.begin(); it != map_of_space_labels.end(); it++){
        auto &label_ = it->first; 
        auto &bunch = it->second;
        std::vector<int> v_clust_labels(bunch.size());
        ClusterizeByVelocity( v_dist, v_method, bunch, v_clust_labels);
        for (size_t i = 0; i < bunch.size(); i++)
            map_of_labels[ std::make_pair(label_, v_clust_labels.at(i))].push_back(bunch.at(i));
    }    
    
    size_t label_counter = 0;
    for (const auto& labl : map_of_labels){
        //auto m_key = labl.first;
        auto &m_parts = labl.second;
        for (const auto& part : m_parts){
            labels_of_particles[part] = label_counter;
        }
        label_counter++;
    }
}

void BunchOfParticles::ClusterizeByVelocityThanBySpace(const double sp_dist, const std::string sp_method, 
                                                       const double v_dist, const std::string v_method){
    if (bunch_of_particles.size() == 0) {
        std::cout << "No particles to clusterize! Exit with doing nothing."<< std::endl;
        return;
    }
    if (bunch_of_particles.size() == 1){
        labels_of_particles = {0};
        return;
    }
    
    std::map< std::pair<int, int>, std::vector<int> > map_of_labels;
    
    std::map <int, std::vector<int>> map_of_vel_labels;    
    ClusterizeByVelocity(v_dist, v_method);    
    for (size_t i = 0; i < labels_of_particles.size(); i++)
        map_of_vel_labels[labels_of_particles.at(i)].push_back(i);
    
    for( auto it = map_of_vel_labels.begin(); it != map_of_vel_labels.end(); it++){
        auto &label_ = it->first;
        auto &bunch = it->second;
        std::vector<int> sp_clust_labels(bunch.size());
        ClusterizeBySpace( sp_dist, sp_method, bunch, sp_clust_labels);
        for (size_t i = 0; i < bunch.size(); i++)
            map_of_labels[ std::make_pair(sp_clust_labels.at(i), label_) ].push_back(bunch.at(i));
    }    
    
    size_t label_counter = 0;
    for (const auto& labl : map_of_labels){
        //auto m_key = labl.first;
        auto &m_parts = labl.second;
        for (const auto& part : m_parts){
            labels_of_particles[part] = label_counter;
        }
        label_counter++;
    }
}

void BunchOfParticles::ClusterizeBySpaceAndVelocityTogether(const double sp_dist, const std::string sp_method, 
                                                            const double v_dist, const std::string v_method){
    if (bunch_of_particles.size() == 0) {
        std::cout << "No particles to clusterize! Exit with doing nothing."<< std::endl;
        return;
    }
    if (bunch_of_particles.size() == 1){
        labels_of_particles = {0};
        return;
    }
    
    ClusterizeByVelocity(v_dist, v_method);
    auto v_clust_labels = labels_of_particles;
    ClusterizeBySpace(sp_dist, sp_method);
    auto sp_clust_labels = labels_of_particles;
    
    std::map< std::pair<int, int>, std::vector<int> > map_of_labels;
    for (size_t i = 0; i < labels_of_particles.size(); i++)
        map_of_labels[ std::make_pair(sp_clust_labels.at(i), v_clust_labels.at(i)) ].push_back(i);
    
    size_t label_counter = 0;
    for (const auto& labl : map_of_labels){
        //auto m_key = labl.first;
        auto &m_parts = labl.second;
        for (const auto& part : m_parts)
            labels_of_particles[part] = label_counter;
        label_counter++;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// General clusterizing function ////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void BunchOfParticles::Clusterize(const double sp_dist, const std::string sp_method, const double v_dist, 
                                  const std::string v_method, const short clust_order){
    switch (clust_order){
        case 0:
            ClusterizeBySpaceThanByVelocity(sp_dist, sp_method, v_dist, v_method);
            break;
        case 1:
            ClusterizeByVelocityThanBySpace(sp_dist, sp_method, v_dist, v_method);
            break;
        case 2:
            ClusterizeBySpaceAndVelocityTogether(sp_dist, sp_method, v_dist, v_method);
            break;
        default:
            std::cout << "Clusterize: Unexpected clust_order value. Allowed values: 0, 1, 2.\nDo nothing."<< std::endl;
            break;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// Get particle labels or 'bunch_of_particles' itself //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const std::vector<int> BunchOfParticles::GetLabelsVector(){
    return labels_of_particles;
}

const std::map<int, std::vector<SingleParticle>> BunchOfParticles::GetLabelsMap(){
    std::map<int, std::vector<SingleParticle>> result;
    for (size_t i = 0; i < labels_of_particles.size(); i++)
        result[labels_of_particles.at(i)].push_back(bunch_of_particles.at(i));
    return result;
}

const std::vector<SingleParticle> BunchOfParticles::GetBunchOfParticles(){
    return bunch_of_particles;
}
