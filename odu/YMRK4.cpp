#include <iostream>
#include <vector>
#include <cmath>

using namespace std;


double pow2 (double x, int n) {
    double res = x;
    for (int i=0; i<n-1; i++) {
        res *= x;
    }
    return res;
}

double f (double x, double u) {
    double func = x + u;
    return func;
}

double calc_new_value (double x, double u0, double h) {
    vector<double> k(4);
    k[0] = f(x, u0);
    k[1] = f(x+h/2, u0+k[0]*h/2);
    k[2] = f(x+h/2, u0+k[1]*h/2);
    k[3] = f(x+h, u0+k[2]*h);
    double u = u0 + h/6*(k[0]+2*k[1]+2*k[2]+k[3]);
    return u;
}

class YMRK4 {
public:
    YMRK4 () {}
    YMRK4 (double begin, double end, double step, double begin_val) { // конструктор для 1 уравнения
        double u0 = begin_val;
        val.push_back(u0);
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
        while (abs(x-(end-begin)) > 1e-6) { // 1/h - кол-во узлов
            double u1, u2, er, u;
            u1 = calc_new_value (x, u0, h);
            u2 = calc_new_value (x, u0, h/2);
            u2 = calc_new_value (x+h/2, u2, h/2);
            er = abs((u2-u1)/(1/pow2(2,p)-1));
            u = u2 - er/pow2(2, p);
            cout << "u1 = " <<u1<<" u2 = "<<u2<<" u = "<<u<<endl;
            if (er <= E) {
                val.push_back(u);
                x +=h;
                nodes.push_back(x);
            } else {
                h = 0.97*pow(E/er, 1./p)*h; // берём на 3% меньше
                cout << "step " << h << endl;
            }
            u0 = u;
        }
    }

    double get_value (double x) { // значени только в узловых точках
        int l = int(x/h);
        return val[l];
    }
    double get_value_any (double x) { // любую точку посчитает
        double y;
        for (int i=0; i<n-1; i++) {
            if (x>=nodes[i] && x<nodes[i+1]) {
                double tau = x-nodes[i];
                y = calc_new_value(nodes[i], val[i], tau);
                break;
            } else if (x >= nodes[n-1]) {
                double tau = x-nodes[n-1];
                y = calc_new_value(nodes[n-1], val[n-1], tau);
            }
        }
        return y;
    }

private:
    vector<double> nodes;
    vector<double> val;
    double h;
    int n;
};


int main (int argc, char* argv[]) {
    YMRK4 m(0., 1., 1e-1, 0);
    double u = m.get_value_any(1);
    cout << "u " << u;
    return 0;
}
