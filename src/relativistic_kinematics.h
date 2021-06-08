#pragma once

#include <cmath>

#include "linear_algebra.h"


std::vector<double> VelocityInNewFrame(const std::vector<double>& v, const std::vector<double>& w){
    //const double c = 299792458;
    const double c = 1.;
    const double w2 = ScalarProduct(w,w); 
    if (ScalarProduct(v,v) > c*c || w2 > c*c){
        std::cerr << "VelocityInNewFrame: One of the velocities is greater than speed of light.\n";
        throw ;
        return {};
    }
    const double vw = ScalarProduct(v,w);
    const double gamma = 1/sqrt(1-w2/c/c); 
    return (1/(1 - vw/c/c)) * ( (1/gamma)*v - w + ((1/c/c)*(gamma/(gamma+1))*vw)*w );
}

std::vector<double> CoordinateInNewFrame(const std::vector<double>& x, const std::vector<double>& v){
    //const double c = 299792458;
    const double c = 1.;
    const double v2 = ScalarProduct(v,v); 
    if (v2 > c*c){
        std::cerr << "CoordinateInNewFrame: Velocity of particle is greater than speed of light.\n";
        throw ;
        return {};
    }
    const double gamma = 1/sqrt(1-v2/c/c);
    const double xv = ScalarProduct(x,v); 
    return x+((gamma-1)*xv/v2)*v;
}
