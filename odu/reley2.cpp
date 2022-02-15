#include <iostream>
#include <vector>
#include <cmath>

#include <fstream>

double f1(double x, double y) {
    return y;
}
double f2(double x, double y) {
    double m = 1000.;
    return m*(1-y*y)*y-x;
}

int main (int argc, char* argv[]) {
    double  h = 1;
    double T = 1000;
    double x1,x2,y1,y2;
    std::vector<double> x(T+1),y(T+1),t(T+1);
    for (int i=0; i<T+1; i++) {
        t[i] = i*h;
    }
    x[0] = 0;
    y[0] = 1e-3;
    for(int i=0; i<int(T); i++) {
        x1 = x[i];
        y1 = y[i];
        x2 = x[i] +h/2*(f1(x[i], y[i])+f1(x1, y1));
        y2 = y[i] +h/2*(f2(x[i], y[i])+f2(x1, y1));
        while (fabs(x2-x1)>1e-6 || fabs(y2-y1)>1e-6) {
            //std::cout << x1 << " "<<x2 << std::endl;
            x1 = x2;
            y1 = y2;
            x2 = x[i] +h/2*(f1(x[i], y[i])+f1(x1, y1));
            y2 = y[i] +h/2*(f2(x[i], y[i])+f2(x1, y1));
            //std::cout << x1 << " "<<x2 << std::endl;
        }
        //std::cout<<"++++++++++++++++++++++++" <<std::endl;
        x[i+1] = x2;
        y[i+1] = y2;

    }
    for(int i=0; i<T+1; i++) {
        std::cout <<t[i]<<" "<< x[i] << std::endl;
    }


    return 0;
}
