#ifndef UTILITIES_PACS1
#define UTILITIES_PACS1


#include<iostream>
#include<functional>
#include<cmath>

struct Parameters{
    unsigned short int dim = 2;
    unsigned int nmax_it = 100;
    unsigned method = 0; //the default one is armijo, o.w. use 1 for exp and 2 for inv
    double sigma = 0.25;
    double mu = 0.2;
    double tol_fun = 1e-6;
    double tol_res = 1e-6;
    double initial_step = 2;
    std::function<double(const std::vector<double> &)> fun;
    std::function<std::vector<double>(const std::vector<double> &)> dfun;
    std::vector<double> initial;
};




// enum step{
//     exp,
//     inv,
//     arm
// };

// template <step lr_choice> 







#endif