#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

double f (double x) { // вычисление подынтегральной функции
    return cos(x)/(2+x*x);
}


double dev2 (double x) { // подсчёт второй производной
    double h = 1e-5;
    return (f(x+h)-2*f(x)+f(x-h))/(h*h);
}


double integr (double h) { // метод трапеций
    double n = 2000/h;
    double I=0;
    for (int i=0; i<n; i++) {
        double x = i*h;
        I+=f((x+x+h)/2)*h;
    }
    return I;
}

int main(int argc, char* argv[]) {
    double h1 = 1e-3;
    double h2 = h1/2;
    double I1 = integr(h1);
    double I2 = integr(h2);
    double I = I2 + (I2-I1)/(2*2-1);
    double E = abs(I1-I2)/(1-0.25);
    cout << I << endl << E << endl << I1;
    return 0;
}
