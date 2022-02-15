#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

double pow2 (double x, int n) {
    double res = x;
    for (int i=0; i<n-1; i++) {
        res *= x;
    }
    return res;
}

struct func {
    double f1,f2,f3,f4;
};

func calc (double x, double u, double y, double v, double h) {
    double m = 0.012277471;
    double et = 1 - m;
    double A = sqrt(pow2(pow2(x+m, 2) + pow2(y, 2), 3));
    double B = sqrt(pow2(pow2(x-et, 2) + pow2(y, 2), 3));
    func s;
    s.f1 = u;
    s.f2 = x + 2*v - et/A*(x+m) - m/B*(x-et);
    s.f3 = v;
    s.f4 = y - 2*u - et*y/A - m*y/B;
    return s;
}


int main (int argc, char* argv[]) {
    std::ofstream file ("/home/senseye3/Study/c_math/odu/out_task_5_NMRK2_10_step.csv");

    double  h = 1e-3;
    double tm = 17.0652165601579625588917206249;
    int T = round(10*tm/h);
    double x1,x2,u1,u2,y1,y2,v1,v2;
    std::vector<double> x(T+1), u(T+1), y(T+1), v(T+1), t(T+1);
    for (int i=0; i<T+1; i++) {
        t[i] = i*h;
    }

    x[0] = 0.994;
    y[0] = 0;
    u[0] = 0;
    v[0] = -2.00158510637908252240537862224;

    for(int i=0; i<T; i++) {
        x1 = x[i];
        u1 = u[i];
        y1 = y[i];
        v1 = v[i];
        func s1 = calc(x[i], u[i], y[i], v[i], h);
        func s2 = calc(x1, u1, y1, v1, h);
        x2 = x[i] + h/2*(s1.f1+s2.f1);
        u2 = u[i] + h/2*(s1.f2+s2.f2);
        y2 = y[i] + h/2*(s1.f3+s2.f3);
        v2 = v[i] + h/2*(s1.f4+s2.f4);

        while (fabs(x2-x1)>1e-6 || fabs(y2-y1)>1e-6) {
            x1 = x2;
            u1 = u2;
            y1 = y2;
            v1 = v2;
            s2 = calc(x1, u1, y1, v1, h);
            x2 = x[i] + h/2*(s1.f1+s2.f1);
            u2 = u[i] + h/2*(s1.f2+s2.f2);
            y2 = y[i] + h/2*(s1.f3+s2.f3);
            v2 = v[i] + h/2*(s1.f4+s2.f4);
        }
        x[i+1] = x2;
        u[i+1] = u2;
        y[i+1] = y2;
        v[i+1] = v2;

    }
    for(int i=0; i<T+1; i++) {
    }
    file << "x,y" << std::endl;
    for (int index = 0; index < t.size(); ++index){
        file << x[index] << "," << y[index] << std::endl;
    }
    file.close();

    return 0;
}
