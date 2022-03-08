#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

struct Area_size {
    double val;
    explicit Area_size (double new_val) {
        val = new_val;
    }
};
struct Step_size {
    double val;
    explicit Step_size (double new_val) {
        val = new_val;
    }
};
struct Courant {
    double val;
    explicit Courant (double new_val) {
        val = new_val;
    }
};
struct Calc_time {
    int val;
    explicit Calc_time (int new_val) {
        val = new_val;
    }
};

double func_initial (double x, int L) {
    return sin(4*M_PI*x/L);
}
double analitic_solution (double x, double t, double L) {
    return sin(4*M_PI*(x-t)/L);
}

class Transfer_eq1 {
public:
    Transfer_eq1 (Area_size new_size, Step_size new_step, Courant new_Co, Calc_time new_time) {
        area_size = new_size.val;
        h = new_step.val;
        NX = int(round(area_size/h))+1;
        Co = new_Co.val;
        tau = Co*h;
        T = new_time.val;
        NT = int(round(T/tau))+1;
        u.resize(NT);
        for (auto& i : u) i.resize(NX);
        //std::cout << u.size() << " " << u[1].size()<< " " << NX;
        nodes.resize(NX);
    }
    void gen_mesh () {
        for (int j=0; j<NX; j++) {
            nodes[j] = j*h;
        }
    }
    void initial_conditions () {
        for (int i=0; i<NX; i++) {
            u[0][i] = func_initial(nodes[i], area_size);
        }
    }
    void lax_ven_solve () {
        for (int i=1; i<NT; i++) {
            for (int j=1; j<NX-1; j++) {
                u[i][j] = (1-Co*Co)*u[i-1][j] + 0.5*(Co*Co-Co)*u[i-1][j+1] + 0.5*(Co*Co+Co)*u[i-1][j-1];
            }
            u[i][0] = (1-Co*Co)*u[i-1][0] + 0.5*(Co*Co-Co)*u[i-1][1] + 0.5*(Co*Co+Co)*u[i-1][NX-1]; // периодическое граничное условие
            u[i][NX-1] = (1-Co*Co)*u[i-1][NX-1] + 0.5*(Co*Co-Co)*u[i-1][0] + 0.5*(Co*Co+Co)*u[i-1][NX-2];
        }
    }
    void solve () {
        for (int i=1; i<NT; i++) {
            for (int j=1; j<NX; j++) {
                u[i][j] = u[i-1][j] - Co*(u[i-1][j] - u[i-1][j-1]); //(u[i-1][j-1] - u[i-1][j])
            }
            u[i][0] = u[i][NX-1]; // периодическое граничное условие
        }
    }

    std::vector<double> get_nodes () {
        return nodes;
    }
    std::vector<std::vector<double>> get_solution () {
        return u;
    }

private:
    double area_size;
    double h;
    int NX;
    double Co;
    double tau;
    int T;
    int NT;
    std::vector<double> nodes;
    std::vector<std::vector<double>> u;
};

int main () {
    std::ofstream file("/home/senseye3/Study/c_math/pde/task1_transfer/problem1_CFL_1.1_5graphs.csv");
    std::ofstream file2("/home/senseye3/Study/c_math/pde/task1_transfer/problem1_CFL_1.1_t_0.csv");
    std::ofstream file3("/home/senseye3/Study/c_math/pde/task1_transfer/problem1_CFL_1.1_t_5.csv");
    std::ofstream file4("/home/senseye3/Study/c_math/pde/task1_transfer/problem1_CFL_1.1_t_10.csv");
    std::ofstream file5("/home/senseye3/Study/c_math/pde/task1_transfer/problem1_CFL_1.1_t_15.csv");
    std::ofstream file6("/home/senseye3/Study/c_math/pde/task1_transfer/problem1_CFL_1.1_t_18.csv");
    Transfer_eq1 t1 = {Area_size(20.), Step_size(0.5), Courant(1.1), Calc_time(18)};
    double tau = 0.5*1.1;
    t1.gen_mesh();
    t1.initial_conditions();
    t1.lax_ven_solve();
    std::vector<double>  x = t1.get_nodes();
    std::vector<std::vector<double>> y = t1.get_solution();
    file << "x(CFL=1.1),y1,y2,y3,y4,y5" << std::endl;
    file2 << "x,y1(time=0_CFL=1.1),y_theory" << std::endl;
    file3 << "x,y1(time=5_CFL=1.1),y_theory" << std::endl;
    file4 << "x,y1(time=10_CFL=1.1),y_theory" << std::endl;
    file5 << "x,y1(time=15_CFL=1.1),y_theory" << std::endl;
    file6 << "x,y1(time=18_CFL=1.1),y_theory" << std::endl;
    for (int index = 0; index < x.size(); ++index){
        file2 << x[index] << "," << y[0][index] << "," << analitic_solution(x[index], 0., 20.) << std::endl;
        file3 << x[index] << "," << y[int(round(5/tau))][index] << "," << analitic_solution(x[index], 5., 20.) << std::endl;
        file4 << x[index] << "," << y[int(round(10/tau))][index] << "," << analitic_solution(x[index], 10., 20.) << std::endl;
        file5 << x[index] << "," << y[int(round(15/tau))][index] << "," << analitic_solution(x[index], 15., 20.) << std::endl;
        file6 << x[index] << "," << y[int(round(18/tau))][index] << "," << analitic_solution(x[index], 18., 20.) << std::endl;
        file << x[index] << "," << y[0][index] << "," << y[int(round(5/tau))][index] << "," << y[int(round(10/tau))][index]
         << "," << y[int(round(15/tau))][index] << "," << y[int(round(18/tau))][index] << std::endl;
    }
    file2.close();
    file3.close();
    file4.close();
    file5.close();
    file6.close();
    file.close();
    return 0;
}
