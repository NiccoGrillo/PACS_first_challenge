#include "solver.cpp"
#include "GetPot"


// manually setting the func and gradient from the main while I fix muparserx
double obj_func(std::vector<double> x){
    return x[0]*x[1] + 4*pow(x[0], 4) + pow(x[1], 2) + 3*x[0];
    //return x[0]*x[0] + x[1]*x[1] - 5;
}
std::vector<double> grad_fun(std::vector<double> x){
    std::vector<double> buf;

    buf.push_back(x[1] + 16*pow(x[0], 3) + 3);
    buf.push_back(x[0] + 2*x[1]);
    return buf;
}






int main(int argc, char **argv){
    

    
    Parameters p;



    GetPot command_line(argc, argv);
    const std::string filename = command_line.follow("parameters_input.dat", 2, "-f", "--file");

    GetPot            gp(filename.c_str());
    

    const std::string section1 = "solver/funcs_and_points/";
    const std::string section2 = "solver/other_params/";

    

    p.dim = gp((section1 + "dim").data(), 2);


    std::string fun_str;
    std::vector<std::string> dfun_str(p.dim, "");
    fun_str = gp((section1 + "fun").data(), "");
    p.initial = std::vector<double>(p.dim, 0.0);




    for(int i = 1; i < p.dim + 1; i++){
        p.initial[i-1] = gp((section1 + "x0_" + std::to_string(i)).data(), 0.0);
        dfun_str[i-1] = gp((section1 + "dfun_" + std::to_string(i)).data(), "");

    }
    p.initial_step = gp((section2 + "iniitial_step").data(), 0.2);
    p.method = gp((section2 + "method").data(), 0);
    p.nmax_it = gp((section2 + "nmax_it").data(), 1e-6);
    p.tol_res = gp((section2 + "tol_res").data(), 1e-6);
    p.sigma = gp((section2 + "sigma").data(), 0.25);
    p.h = gp((section2 + "h").data(), 0.1);
    p.mu = gp((section2 + "mu").data(), 0.2);


    std::cout<<"MAIN FILE: "<< p.method<< std::endl;

    //for now initialize the functions like this:
    p.fun = obj_func;
    p.dfun = grad_fun;


    //Initialize the solver and compute the solution

    Solver my_solver(p);

    my_solver.solve();



    return 0;
    



    


}