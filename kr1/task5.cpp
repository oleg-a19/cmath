#include <iostream>
#include <cmath>

using namespace std;

double pow2 (double x) {
    if (x >= 0) return pow(x,1./3.);
    else return -pow(-x,1./3.);
}

int main(int argc, char* argv[]) {
    double y0=0.91;
    double y02 = 0.9;
    double d = 1e-2;
    double y2 = pow2(-(-(3+d)+(2+d)*y02*y02)/(1+d));
    double y = sqrt(exp(y0)/3);
    while (abs(y-y0) > 0.0001) {
        y0 = y;
        y = sqrt(exp(y0)/3);
    }
    cout << fixed;
    cout.precision(10);
    cout << "y = " << y ;
    return 0;
}
