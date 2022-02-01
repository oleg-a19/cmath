#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;

double f (double x) { // вычисление подынтегральной функции
    return x*x*x+5*x*x+10*x+4.5;
}

int main(int argc, char* argv[]) {
    double x1 = 0.5*(-3*sqrt(15)/5-1);
    double x2 = -0.5;
    double x3 = 0.5*(3*sqrt(15)/5-1);
    double w1 = 5./18;
    double w2 = 4./9;
    double w3 = w1;
    double I = 3*(w1*f(x1)+w2*f(x2)+w3*f(x3));
    cout << I;

    return 0;
}
