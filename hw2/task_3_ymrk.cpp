#include <iostream>
#include <vector>
#include <cmath>
#include<iomanip>
#include <fstream>

using namespace std;

// замечания: конструктор лёгкий(все параметры передать в метод класса), вектор структур(х, у, у'), методы доступа к полям, не писать using

double pow2 (double x, int n) {
    double res = x;
    for (int i=0; i<n-1; i++) {
        res *= x;
    }
    return res;
}

struct function {
    double f1, f2;
};

function f (double x, double y, double z) {
    function s;
    s.f1 = z;
    s.f2 = (2-6*x+2*pow2(x,3) + (x*x-3)*exp(x)*sin(x)*(1+cos(x))+
        cos(x)*(exp(x)+(x*x-1)+pow2(x, 4)-3*x*x)) - z*(x*x-3) - y*(x*x-3)*cos(x);
    return s;
}

vector<double> calc_new_value (double x, double y, double z, double h) {
    vector<double> k1(4);
    vector<double> U(2);
    k1[0] = f(x, y, z).f1;
    k1[1] = f(x+h/2, y+k1[0]*h/2, z+k1[0]*h/2).f1;
    k1[2] = f(x+h/2, y+k1[1]*h/2, z+k1[1]*h/2).f1;
    k1[3] = f(x+h, y+k1[2]*h, z+k1[2]*h).f1;
    U[0] = y + h/6*(k1[0]+2*k1[1]+2*k1[2]+k1[3]);

    vector<double> k2(4);
    k2[0] = f(x, y, z).f2;
    k2[1] = f(x+h/2, y+k2[0]*h/2, z+k2[0]*h/2).f2;
    k2[2] = f(x+h/2, y+k2[1]*h/2, z+k2[1]*h/2).f2;
    k2[3] = f(x+h, y+k2[2]*h , z+k2[2]*h).f2;
    U[1] = z + h/6*(k2[0]+2*k2[1]+2*k2[2]+k2[3]);
    return U;
}

class YMRK4 {
public:
    YMRK4 () {}
    YMRK4 (double begin, double end, double step, double begin_val, double begin_val2) { // конструктор для 2 уравнений
        double y = begin_val;
        double z = begin_val2;
        solution.push_back(y);
        solution_derivative.push_back(z);
        int p=4; // порядок аппроксимации
        h = step;
        n = int((end-begin)/h)+1;
        double E = 1e-3;
        /*for (int i=0; i<n; i++) {
            nodes.push_back(i*h);
            cout << nodes[i] << endl;
        }*/
        double x = 0;
        nodes.push_back(x);
        while ((end-begin)-x > 1e-6) { // 1/h - кол-во узлов
            vector<double> u1, u2, u(2);
            double er_y, er_z;
            u1 = calc_new_value (x, y, z, h);
            u2 = calc_new_value (x, y, z, h/2);
            u2 = calc_new_value (x+h/2, u2[0], u2[1], h/2);
            er_y = abs((u2[0]-u1[0])/(1/pow2(2,p)-1));
            er_z = abs((u2[1]-u1[1])/(1/pow2(2,p)-1));
            u[0] = u2[0] - er_y/pow2(2, p);
            u[1] = u2[1] - er_z/pow2(2, p);
            //cout << "er_y-E =  " <<er_y-E<<" er_z-E "<<er_z-E<<endl;
            if ((er_y-E)<0 && (er_z-E)<0) {
                //cout << h << endl;
                solution.push_back(u[0]);
                solution_derivative.push_back(u[1]);
                x +=h;
                nodes.push_back(x);
                //cout << abs(x-(end-begin)) << endl;
            } else {
                if (er_y <= er_z) h = 0.97*pow(E/er_z, 1./p)*h; // берём на 3% меньше
                else h = 0.97*pow(E/er_y, 1./p)*h; // берём на 3% меньше
                //cout << "er_y-E =  " <<er_y-E<<" er_z-E = "<<er_z-E<<endl;
                //cout << "step " << h << endl;
            }
            y = u[0];
            z = u[1];
        }
    }
    double get_value (double x) { // значени только в узловых точках
        int l = int(x/h);
        return solution[l];
    }
    double get_value_any (double x) { // любую точку посчитает
        double y;
        for (int i=0; i<n-1; i++) {
            if (x>=nodes[i] && x<nodes[i+1]) {
                double tau = x-nodes[i];
                y = calc_new_value(nodes[i], solution[i], solution_derivative[i], tau)[0];
                break;
            } else if (x >= nodes[n-1]) {
                double tau = x-nodes[n-1];
                y = calc_new_value(nodes[n-1], solution[n-1], solution_derivative[n-1], tau)[0];
            }
        }
        return y;
    }

    vector<double> nodes;
    vector<double> solution, solution_derivative;
    double h;
    int n;
};


int main (int argc, char* argv[]) {
    std::ofstream file ("/home/senseye3/Study/c_math/odu/out_task_3_ymrk.csv");

    YMRK4 m(0., M_PI, 1e-3, 0, 0.50979559);
    for (int i=0; i<m.nodes.size(); i++) {
        //cout <<"x = "<<m.nodes[i]<<" y = "<< m.solution[i]<<endl;
    }
    file << "x,y" << std::endl;
    for (int index = 0; index < m.solution.size(); ++index){
        file << m.nodes[index] << "," << m.solution[index] << std::endl;
    }
    double u = m.get_value_any(M_PI);
    cout << fixed;
    cout << "u " << setprecision(5) << u;
    file.close();
    return 0;
}
