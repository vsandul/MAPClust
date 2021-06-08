#pragma once

#include <vector>

#include "linear_algebra.h"

const int PROTON_PDG = 2212, NEUTERON_PDG = 2112, LAMBDA_PDG = 3122;

const int ClusterPDGCode( const std::vector<int>& partons){
    int nN = 0, nP = 0, nL = 0;
    if (partons.size() == 0) return 0;
    if (partons.size() == 1) return partons.at(0);
    for (size_t i = 0; i < partons.size(); i++){
        if (partons.at(i) == PROTON_PDG) {
            nP++;
            continue;
        }
        if (partons.at(i) == NEUTERON_PDG) {
            nN++;
            continue;
        }
        if (partons.at(i) == LAMBDA_PDG) {
            nL++;
            continue;
        }
    }
    
    int A = nN + nP + nL;
    if (A == 0) return 0;
    if (A == 1) {
        if (nN == 1) return NEUTERON_PDG;
        if (nP == 1) return PROTON_PDG;
        if (nL == 1) return LAMBDA_PDG;
    }
    
    return A*10 + nP*1E4 + nL*1E7 + 1E9; // ±10LZZZAAAI    
}

//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------

template<typename T>
std::vector<T> VectorTotal(const vector<vector<T>>& val){
    vector<T> res(val.size(), 0);
    for (const auto& i:val)
        res = res + i;
    return res;
}

template<typename T>
T ScalarTotal(const vector<T>& val){
    T res = 0;
    for (const auto& i:val)
        res += i;
    return 1.0*res;
}

template<typename T>
std::vector<double> VectorAverage(const vector<vector<T>>& val){
    return (1./val.size())*VectorTotal(val);
}

template<typename T>
double ScalarAverage(const vector<T>& val){
    return (1.0/val.size())*ScalarTotal(val);
}

