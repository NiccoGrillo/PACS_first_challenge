#include "solver.hpp"

Solver::Solver(const Parameters &param_): param(param_){
    solution.resize(param.dim, 0.0);
};



void Solver::solve(){
        std::vector<double> x_k = param.initial;
    double f_eval_k = param.fun(x_k);
    std::vector<double> grad_k(x_k);

    double diff_x_k = 0.0;
    double diff_fun_k = 0.0;

    double buffer_x = 0.0, buffer_fun= param.fun(x_k); //two buffers used later to update the diff_x_k and diff_fun_k values

    double alpha_k = param.initial_step;

    bool converged = false;
        for (unsigned m_iter = 1; !converged; ++m_iter)
        { //evaluate function and gradiant at point x_k
            grad_k = param.dfun(x_k);

            //std::cout<<"iteration: "<< m_iter<<";"<< x_k[0]<<";"<<x_k[1]<< std::endl;

            alpha_k = search_alpha(param.method, x_k, grad_k, buffer_fun, m_iter);//computes alpha_k with expressed method, default: armijo


            //now we update x_k and compute diff_x_k

            diff_x_k = 0.0;
            for (unsigned kk = 0; kk < param.dim; ++kk){
                buffer_x = x_k[kk]; 
                x_k[kk] = x_k[kk] - alpha_k*grad_k[kk];

                diff_x_k += fabs(buffer_x - x_k[kk]);

            }

            f_eval_k = param.fun(x_k);




            //now we update diff_fun_k
            diff_fun_k = fabs(f_eval_k - buffer_fun);
            buffer_fun = f_eval_k;



            //and finally we check convergence
            if (diff_fun_k < param.tol_fun){
                std::cout<<"reached convergence through FUNCTION RESIDUAL, after: " << m_iter << " iterations"<<std::endl;
                converged = true;
            }
            if (diff_x_k <param.tol_res){
                std::cout<<"reached convergence through VALUE RESIDUAL, after: " << m_iter << " iterations"<<std::endl;
                converged = true;
            }
            if (m_iter>= param.nmax_it){
                std::cout<<"reached MAXIMUM number of ITERATIONS: " << m_iter <<std::endl;
                converged = true;
            }

        }
        solution = x_k;
        print_solution();

}


double Solver::search_alpha(const unsigned short method,
            const std::vector<double> &x,
            const std::vector<double> &grad_k,
            const double &val,
            const unsigned k){
    double alpha_k = param.initial_step;


    if (method == 1){//exponential method
        alpha_k = param.initial_step*exp(-1*param.mu * k);
        //std::cout<<"ci arrivo qui?    alpha:"<<param.initial_step << ";"<< alpha_k<< ";"<<param.mu << ";"<< k<<std::endl;
    }
    else if (method == 2){ //inverse method
        alpha_k = param.initial_step * (1/(1+param.mu * k));
    }
    else {  //amijo, which is the default one
        std::vector<double> diff_vec(param.dim, 0);



        for (unsigned kk = 0; kk < param.dim; kk++){
            diff_vec[kk] = x[kk] - alpha_k*grad_k[kk];
        }
        //bool check = true;
        while(val - param.fun(diff_vec) < param.sigma*alpha_k*compute_norm2(grad_k)){
            alpha_k = alpha_k/2.0;
            for (unsigned kk = 0; kk < param.dim; kk++){
                diff_vec[kk] = x[kk] - alpha_k*grad_k[kk];
            }
        }
        }
    return alpha_k;
}
    





double Solver::compute_norm2(const std::vector<double> &vec){
    double sum = 0.0;
    for (auto i:vec){
        sum += i*i;
    }
    return sum;

}