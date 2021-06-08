#include "dsu.h"

//vector realization
int dsu_get (int v, std::vector<int>& p) {
	return (v == p[v]) ? v : (p[v] = dsu_get (p[v], p));
}

void dsu_unite (int a, int b, std::vector<int>& p) {
	a = dsu_get (a, p);
	b = dsu_get (b, p);
	if (rand() & 1)
		std::swap (a, b);
	if (a != b)
		p[a] = b;
}


//array realization
int dsu_get (int v, int* p) {
	return (v == p[v]) ? v : (p[v] = dsu_get (p[v], p));
}

void dsu_unite (int a, int b, int* p) {
	a = dsu_get (a, p);
	b = dsu_get (b, p);
	if (rand() & 1)
		std::swap (a, b);
	if (a != b)
		p[a] = b;
}
