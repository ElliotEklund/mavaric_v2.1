#ifndef lyp_funcs_hpp
#define lyp_funcs_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <cmath>
#include <tgmath.h>

template <typename T>
T p_norm(double p, const vector<T> & v){
    
    int n = v.size();
    T norm(0);
    
    for (int i=0; i<n; i++) {
        norm +=  pow(fabs(v(i)),p);
    }
    
    return norm;
}

#endif
