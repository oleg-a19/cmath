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

double func_initial_1 (double x) {
    return sin(x)+cos(x);
}
double func_initial_2 (double x) {
    return sin(x)-cos(x);
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
        u1.resize(NT);
        u2.resize(NT);
        for (auto& i : u1) i.resize(NX);
        for (auto& i : u2) i.resize(NX);
        //std::cout << u1.size() << " " << u1[1].size()<< " " << NX;
        nodes.resize(NX);
    }
    void gen_mesh () {
        for (int j=0; j<NX; j++) {
            nodes[j] = j*h;
        }
    }
    void boundary_conditions () {
        for (int i=1; i<NT; i++) {
            u1[i][0] = 0;
            u1[i][NX-1] = 0;
            u2[i][0] = 0;
            u2[i][NX-1] = 0;
        }
    }
    void initial_conditions () {
        for (int i=0; i<NX; i++) {
            u1[0][i] = func_initial_1(nodes[i]);
            u2[0][i] = func_initial_2(nodes[i]);
        }

    }
    void solve () {
        for (int i=1; i<NT; i++) {
            if (i >= NX) break;
            for (int j=1; j<NX-i; j++) {
                u1[i][j] = u1[i-1][j] - Co*(u1[i-1][j] - u1[i-1][j-1]);
                u2[i][j] = u2[i-1][j] + Co*(u2[i-1][j] - u2[i-1][j-1]);

            }
        }
    }

    std::vector<double> get_nodes () {
        return nodes;
    }
    std::vector<std::vector<double>> get_solution_1 () {
        /*std::vector<std::vector<double>> v1(NT);
        for (auto& i : v1) i.resize(NX);
        for (int i=0; i < NT; i++) {
            for (int j=0; j<NX; j++) {
                v1[i][j] = 0.5*(u1[i][j]+u2[i][j]);
            }
        }*/
        return u1;
    }
    std::vector<std::vector<double>> get_solution_2 () {
        /*std::vector<std::vector<double>> v2(NT);
        for (auto& i : v2) i.resize(NX);
        for (int i=0; i < NT; i++) {
            for (int j=0; j<NX; j++) {
                v2[i][j] = 0.5*(u1[i][j]+u2[i][j]);
            }
        }*/
        return u2;
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
    std::vector<std::vector<double>> u1;
    std::vector<std::vector<double>> u2;
};

int main () {
    std::ofstream file("/home/senseye3/Study/c_math/pde/sr/t1.csv");
    Transfer_eq1 t1 = {Area_size(M_PI), Step_size(0.05), Courant(0.6), Calc_time(10)};
    t1.gen_mesh();
    t1.boundary_conditions();
    t1.initial_conditions();
    t1.solve();
    std::vector<double>  x = t1.get_nodes();
    std::vector<std::vector<double>> y1 = t1.get_solution_1();
    std::vector<std::vector<double>> y2 = t1.get_solution_2();
    file << "x,y1(time=0),y2(time=0.5),y3(time=1),y4(time=1.5)" << std::endl;
    for (int index = 0; index < x.size(); ++index){
        //std::cout << x[index] << " zzz " << y[0][index] << std::endl;
        file << x[index] << "," << y1[0][index] << "," << y1[int(round(0.5/0.03))][index] << "," << y1[int(round(1/0.03))][index]
         << "," << y1[int(round(1.5/0.03))][index] << std::endl;
    }
    file.close();
    return 0;
}
