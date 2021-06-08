#include "dbscan.h"


std::vector<int> RangeQuery(size_t num_of_part, double* distances,
		const int& partnum, const double& eps ){
	std::vector<int> neighbors;
	for (size_t i = 0; i < num_of_part; i++){
		if (i == partnum || distances[i] >= eps) continue;
		else neighbors.push_back(i);
	}
	return neighbors;
}

std::vector<int> DBSCAN( double** distmat, size_t num_of_part, double eps, short minPts){
	size_t C = -1; // cluster counter
    int *labels = new int[num_of_part];
    for (int i = 0; i < num_of_part; i++)
        labels[i] = UNDEFINED;
	for (int i = 0; i < num_of_part; i++){
		if (labels[i] != UNDEFINED) continue;
		std::vector<int> neighbors = RangeQuery(num_of_part, distmat[i], i, eps);
		if (neighbors.size() < minPts){
			labels[i] = NOISE;
			continue;
		}
		labels[i] = ++C;
		std::vector<int> seedset(neighbors.begin(), neighbors.end());
		for (size_t s = 0; s < seedset.size(); s++){
            int p = seedset.at(s);
			if (labels[p] == NOISE)	labels[p] = C;
			if (labels[p] != UNDEFINED)	continue;
			labels[p] = C;
			neighbors = RangeQuery(num_of_part, distmat[p], p, eps);
			if (neighbors.size() >= minPts)
				seedset.insert(seedset.end(), neighbors.begin(), neighbors.end());            
		}
	}
	
	int noise_counter = 0;
	for (int i = 0; i < num_of_part; i++){
        if (labels[i] == -2)
            labels[i] -= ++noise_counter;
    }        
	return std::vector<int> (labels, labels+num_of_part);
}
