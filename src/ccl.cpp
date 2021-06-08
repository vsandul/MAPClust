#include "ccl.h"

//array realization
std::vector<int> CCLClust(double** distmat, int npoints, double rmax){
	int* labels = new int[npoints];
	for (int i = 0; i < npoints; i++){
		labels[i] = i;
	}

	for (int i = 0; i < npoints; i++){
		for (int j = i+1; j < npoints; j++){
			if(distmat[i][j] < rmax && dsu_get (i,  labels) != dsu_get (j, labels) ){
				dsu_unite(i, j, labels);
			}
		}
	}

	std::map<int, std::vector<int>> lab_map;
	for (size_t i = 0; i < npoints; i++)
		lab_map[dsu_get(labels[i], labels)].push_back(i);
	int lab_counter = 0;
	for (const auto& l:lab_map){
		for (const auto& node:l.second)
			labels[node] = lab_counter;
		lab_counter++;
	}
	return std::vector<int> (labels, labels+npoints);
}
