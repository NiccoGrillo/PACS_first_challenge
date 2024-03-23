#ifndef SOLVER_PACS1
#define SOLVER_PACS1

#include "Utilities.hpp"




class Solver{
    protected:
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

        virtual double search_alpha(
            const std::vector<double> &x,
            const std::vector<double> &grad_k,
            const double &val,
            const unsigned k) = 0;
        
        double compute_norm2(const std::vector<double> &vec);

        

};

class Solver_inv: public Solver{
    public:
        Solver_inv(const Parameters &param_): Solver(param_){};

        double search_alpha(
        const std::vector<double> &x,
        const std::vector<double> &grad_k,
        const double &val,
        const unsigned k) override;
};

class Solver_exp: public Solver{
    public:
        Solver_exp(const Parameters &param_): Solver(param_){};


        double search_alpha(
        const std::vector<double> &x,
        const std::vector<double> &grad_k,
        const double &val,
        const unsigned k) override;
};

class Solver_arm: public Solver{
    public:
        Solver_arm(const Parameters &param_): Solver(param_){};


        double search_alpha(
        const std::vector<double> &x,
        const std::vector<double> &grad_k,
        const double &val,
        const unsigned k) override;
};


#endif