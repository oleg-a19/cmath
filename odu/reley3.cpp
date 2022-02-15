#include <iostream>
#include <vector>
#include <cmath>

#include <fstream>

double f1(double x, double y, double xi, double yi, double h) {
    //double m = 1000.;
    double f = xi + h/2*(yi+y);
    //double f = m*((1-yi*yi)*yi + (1-y*y)*y) - 2/h*(y-yi) - xi;
    return f;
}
double f2(double x, double y, double xi, double yi, double h) {
    double m = 1000.;
    double f = yi + h/2*(m*((1-yi*yi)*yi+(1-y*y)*y)-xi-x);
    //double f = sqrt(1-((2/h*(y-yi)+(xi+x))/m - yi*(1-yi*yi))/y);
    //double f = ((2/h*(y-yi)+(xi+x))/m - yi*(1-yi*yi))/(1-y*y);
    //double f = 2/h*(x-xi)-yi;
    return f;
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
        x2 = f1(x1, y1, x[i], y[i], h);
        y2 = f2(x1, y1, x[i], y[i], h);

        //std::cout << fabs(x2-x1) <<std::endl;// << " " << x2 << " " << y1 << " " << y2 << std::endl;
        int j=0;
        while (fabs(x2-x1)>1e-6 || fabs(y2-y1)>1e-6) {
            //std::cout<<"+++"<<std::endl;
            //j++;
            x1 = x2;
            y1 = y2;
            x2 = f1(x1, y1, x[i], y[i], h);
            y2 = f2(x1, y1, x[i], y[i], h);
            std::cout << x1 << " " << x2 << " " << y1 << " " << y2 << std::endl;
            //std::cout << x1-x2 <<std::endl;
            if (j==10) break;
        }
        //std::cout<<"+++"<<std::endl;
        if (j==10) break;
        x[i+1] = x2;
        y[i+1] = y2;

    }
    for(int i=0; i<T+1; i++) {
        //std::cout <<t[i]<<" "<< x[i] << std::endl;
    }


    return 0;
}
