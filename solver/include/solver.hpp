#ifndef SOLVER_PACS1
#define SOLVER_PACS1

#include "Utilities.hpp"




class Solver{
    private:
        Parameters param;
        std::vector<double> solution;
        
        void print_solution(){
            std::cout<<"converged to the vector: { ";
            for (auto i:solution){
                std::cout<<" "<<i<<" ";
            }
            std::cout<<"}"<<std::endl;
        }
    public:
        Solver(const Parameters &param_);

        void solve();

        double search_alpha(const unsigned short method,
            const std::vector<double> &x,
            const std::vector<double> &grad_k,
            const double &val,
            const unsigned k);
        
        double compute_norm2(const std::vector<double> &vec);

        

};


#endif