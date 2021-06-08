#pragma once

#include <vector>


enum DBSCANMarks{
	UNDEFINED = -1,
	NOISE = -2,
};

std::vector<int> RangeQuery(size_t num_of_part, double* distances,
		const int& partnum, const double& eps );
std::vector<int> DBSCAN( double** distmat, size_t num_of_part, double eps, short minPts);

