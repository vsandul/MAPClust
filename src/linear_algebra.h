#pragma once

#include <vector>
#include <iostream>

template <typename T>
T ScalarProduct(const std::vector<T>& x, const std::vector<T>& y){
	if (x.size() != y.size()){
        std::cerr << "ScalarProduct: Vectors have different sizes\n";
        throw ;
        return (double)NULL;
    }
	if (x.size() == 0 || y.size() == 0){
        std::cerr << "ScalarProduct: Vectors have zero length\n";
		throw ;
        return (double)NULL;
	}
    T result = 0;
    for (size_t i = 0; i < x.size(); i++)
        result += x.at(i)*y.at(i);
    return result;        
}

template <typename T>
T AbsValue(const std::vector<T>& x){
	if (x.size() == 0 ){
        std::cerr << "AbsValue: Vector has zero length\n";
		throw ;
        return (double)NULL;
	}
    T result = 0;
    for (size_t i = 0; i < x.size(); i++)
        result += x.at(i)*x.at(i);
    return result;        
}

template<typename T1, typename T2>
std::vector<T2> operator * (const T1& a, const std::vector<T2>& x){
	std::vector<T2> result(x.size());
	if (x.size() != 0) {
		for (size_t i = 0; i < x.size(); i++){
			result[i] = x.at(i)*a;
		}
	}
	return result;
}

template<typename T1, typename T2>
std::vector<T1> operator * (const std::vector<T1>& x, const T2& a){
	std::vector<T1> result(x.size());
	if (x.size() != 0) {
		for (size_t i = 0; i < x.size(); i++){
			result[i] = x.at(i)*a;
		}
	}
	return result;
}

template<typename T>
std::vector<T> operator - (const std::vector<T>& x, const std::vector<T>& y){
    if (x.size() != y.size()){
        std::cerr << "operator - : Vectors have different sizes\n";
        throw ;
        return {};
    }
    std::vector<T> w(x.size());
    for(size_t i = 0; i < x.size(); ++i)
        w[i] = x.at(i) - y.at(i);
    return w;
}

template<typename T>
std::vector<T> operator + (const std::vector<T>& x, const std::vector<T>& y){
    if (x.size() != y.size()){
        std::cerr << "operator - : Vectors have different sizes\n";
        throw ;
        return {};
    }
    std::vector<T> w(x.size());
    for(size_t i = 0; i < x.size(); ++i)
        w[i] = x.at(i) + y.at(i);
    return w;
}
