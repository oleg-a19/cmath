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
struct Coef_transfer {
    std::vector<double> val;
    explicit Coef_transfer (std::vector<double> new_val) {
        val = new_val;
    }
};
struct Right_part {
    std::vector<double> val;
    explicit Right_part (std::vector<double> new_val) {
        val = new_val;
    }
};

double func_initial_1 (double x) {
    return sin(x)+cos(x);
}
double func_initial_2 (double x) {
    return sin(x)-cos(x);
}
double solution_1_theory_b (double x, double t) {
    return (sin(t)+cos(t))*sin(x) + t;
}
double solution_2_theory_b (double x, double t) {
    return (cos(t)-sin(t))*cos(x);
}

class Transfer_eq1 {
public:
    Transfer_eq1 (Area_size new_size, Step_size new_step, Courant new_Co, Calc_time new_time, Coef_transfer new_coef, Right_part new_part) {
        area_size = new_size.val;
        h = new_step.val;
        Co = new_Co.val;
        T = new_time.val;
        a = new_coef.val;
        f = new_part.val;
        tau = Co*h;
        NX = int(round(area_size/h))+1;
        NT = int(round(T/tau))+1;
        u1.resize(NT);
        u2.resize(NT);
        for (auto& i : u1) i.resize(NX);
        for (auto& i : u2) i.resize(NX);
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
            for (int j=i; j<NX-i; j++) {
                u1[i][j] = f[0]*tau + u1[i-1][j] - Co/2*((a[0]+fabs(a[0]))*(u1[i-1][j]-u1[i-1][j-1]) + (a[0]-fabs(a[0]))*(u1[i-1][j+1]-u1[i-1][j]));
                u2[i][j] = f[1]*tau + u2[i-1][j] - Co/2*((a[1]+fabs(a[1]))*(u2[i-1][j]-u2[i-1][j-1]) + (a[1]-fabs(a[1]))*(u2[i-1][j+1]-u2[i-1][j]));
                //u1[i][j] =tau*f[0] + u1[i-1][j] - Co*(u1[i-1][j] - u1[i-1][j-1]); // правый уголок
                //u2[i][j] =tau*f[1] + u2[i-1][j] + Co*(u2[i-1][j+1] - u2[i-1][j]); // левый уголок
            }
        }
    }

    std::vector<double> get_nodes () {
        return nodes;
    }
    std::vector<std::vector<double>> get_solution_1 () {
        std::vector<std::vector<double>> v1(NT);
        for (auto& i : v1) i.resize(NX);
        for (int i=0; i < NT; i++) {
            for (int j=0; j<NX; j++) {
                v1[i][j] = 0.5*(u1[i][j]+u2[i][j]);
            }
        }
        return v1;
    }
    std::vector<std::vector<double>> get_solution_2 () {
        std::vector<std::vector<double>> v2(NT);
        for (auto& i : v2) i.resize(NX);
        for (int i=0; i < NT; i++) {
            for (int j=0; j<NX; j++) {
                v2[i][j] = 0.5*(u1[i][j]-u2[i][j]);
            }
        }
        return v2;
    }

private:
    double area_size;
    double h;
    int NX;
    double Co;
    double tau;
    int T;
    int NT;
    std::vector<double> a;
    std::vector<double> f;
    std::vector<double> nodes;
    std::vector<std::vector<double>> u1;
    std::vector<std::vector<double>> u2;
};

int main () {

    Transfer_eq1 t1 = {Area_size(M_PI), Step_size(0.05), Courant(0.6), Calc_time(10), Coef_transfer({1., -1.}), Right_part({1., 1.})};
    t1.gen_mesh();
    t1.boundary_conditions();
    t1.initial_conditions();
    t1.solve();
    std::vector<double>  x = t1.get_nodes();
    std::vector<std::vector<double>> y1 = t1.get_solution_1();
    std::vector<std::vector<double>> y2 = t1.get_solution_2();

    std::ofstream file("/home/senseye3/Study/c_math/pde/sr/u2_b_time_0.csv");
    std::ofstream file2("/home/senseye3/Study/c_math/pde/sr/u2_b_time_0.5.csv");
    std::ofstream file3("/home/senseye3/Study/c_math/pde/sr/u2_b_time_0.9.csv");
    file << "x,y2(time=0)_b, y2_theory_b" << std::endl;
    file2 << "x,y2(time=0.5)_b, y2_theory_b" << std::endl;
    file3 << "x,y2(time=0.9)_b, y2_theory_b" << std::endl;
    std::ofstream file4("/home/senseye3/Study/c_math/pde/sr/u2_b_time_0.0-0.9.csv");
    file4 << "x,y2(b)_(time=0)_b,y2(b)_(time=0.5),y2(b)_(time=0.9)" << std::endl;

    for (int index = 0; index < x.size(); ++index){
        file << x[index] << "," << y2[0][index] << "," << solution_2_theory_b(x[index], 0.) << std::endl;
        file2 << x[index] << "," << y2[int(round(0.5/0.03))][index] << "," << solution_2_theory_b(x[index], 0.5) << std::endl;
        file3 << x[index] << "," << y2[int(round(0.9/0.03))][index] << "," << solution_2_theory_b(x[index], 0.9) << std::endl;
        file4 << x[index] << "," << y2[int(round(0./0.03))][index] << "," << y2[int(round(0.5/0.03))][index] << "," << y2[int(round(0.9/0.03))][index] << std::endl;
    }

    file2.close();
    file3.close();
    file.close();
    file4.close();
    return 0;
}
