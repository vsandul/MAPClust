#pragma once

#include <iostream>
#include <vector>
#include <map>
#include "float.h"
#include "dsu.h"

const double M = DBL_MAX;

//array realization
std::vector<int> CCLClust(double** distmat, int npoints, double rmax);
