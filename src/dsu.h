#pragma once

#include <vector>
#include <cstdlib>

//vector realization
int dsu_get (int v, std::vector<int>& p);
void dsu_unite (int a, int b, std::vector<int>& p);

//array realization
int dsu_get (int v, int* p);
void dsu_unite (int a, int b, int* p);
