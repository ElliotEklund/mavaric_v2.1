#ifndef simpson_hpp
#define simpson_hpp

#include <boost/numeric/ublas/vector.hpp>

/* Compute numerical integral using Simpson's method for evenly spaced,
 one-dimensional integration.
 
 delta: width of grid
 n: number of grid points. must be odd
 f: vector of samples evaluated over grid */
template <typename T>
T simp(const double delta, const int n, const vector<T> &f){

    int m = (n-1)/2;
    double sum_odd = 0;
    double sum_even = 0;
    double sum_tails = f(0) + f(n-1);
        
    //Odd sum
    for(int i=1; i<=m; i++){
        sum_odd += f(2*(i-1) + 1);
    }
    
    //Even sum
    for(int i=1; i<=m-1; i++){
        sum_even += f(2*i);
    }
    
    return delta*(sum_tails + 4.0*sum_odd + 2.0*sum_even)/3.0;
}

#endif
