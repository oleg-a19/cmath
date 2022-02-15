#include <iostream>
#include <vector>
#include <cmath>
#include "run_through_3d.h"
#include <fstream>

double pow2 (double x, int n) {
    double res = x;
    for (int i=0; i<n-1; i++) {
        res *= x;
    }
    return res;
}

struct function {
    std::vector<double> x, a, b, c, d;
};

function calc (double h, double first_node, double last_node) {
    function f;
    int n = round((last_node-first_node)/h)+1;
    f.x.resize(n);
    f.a.resize(n);
    f.b.resize(n);
    f.c.resize(n);
    f.d.resize(n);
    f.b[0] = 1;
    f.b[n-1] = 1;
    f.c[0] = 0;
    f.c[n-1] = 0;
    f.a[n-1] = 0;
    f.a[0] = 0;
    f.d[0] = 0;
    f.d[n-1] = M_PI*M_PI;
    f.x[0] = first_node;
    f.x[n-1] = last_node;
    for (int i=1; i<n-1; ++i) {
        f.x[i] = i*h;
        //f.a[i] = 1 + h*(f.x[i]*f.x[i]-3);
        f.a[i] = 1 + h/2*(f.x[i]*f.x[i]-3);
        //f.b[i] = -2 - h*(f.x[i]*f.x[i]-3) + h*h*(f.x[i]*f.x[i]-3)*cos(f.x[i]);
        f.b[i] = -2 + h*h*(f.x[i]*f.x[i]-3)*cos(f.x[i]);
        //f.c[i] = 1;
        f.c[i] = 1 - h/2*(f.x[i]*f.x[i]-3);
        if ((fabs(f.b[i])-fabs(f.a[i])-fabs(f.c[i]))<0) std::cout <<"+++"<< endl;
        f.d[i] = h*h*(2-6*f.x[i]+2*pow2(f.x[i],3) + (f.x[i]*f.x[i]-3)*exp(f.x[i])*sin(f.x[i])*(1+cos(f.x[i]))+
            cos(f.x[i])*(exp(f.x[i])+(f.x[i]*f.x[i]-1)+pow2(f.x[i], 4)-3*f.x[i]*f.x[i]));
    }
    return f;
}

int main (int argc, char* argv[]) {
    std::ofstream file ("/home/senseye3/Study/c_math/odu/out_task_3.csv");
    double h = 1e-3;
    double x0 = 0;
    double x_n = M_PI;
    function f = calc(h, x0, x_n);
    run_through solution(f.a, f.b, f.c, f.d);
    std::vector<double> y = solution.get_solution();

    file << "x,y" << std::endl;
    for (int index = 0; index < f.x.size(); ++index){
        file << f.x[index] << "," << y[index] << std::endl;
    }
    file.close();
    return 0;
}
