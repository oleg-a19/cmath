#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

double f (double x) { // вычисление подынтегральной функции
    return sin(100.*x)*exp(-x*x)*cos(2.*x);
}


double dev2 (double x) { // подсчёт второй производной
    double h = 1e-5;
    return (f(x+h)-2*f(x)+f(x-h))/(h*h);
}


double integr (double h) { // метод трапеций
    double n = 3/h +1;
    double I=0;
    for (int i=0; i<n-1; i++) {
        double x = i*h;
        I+=(f(x)+f(x+h))/2*h;
    }
    return I;
}

int main(int argc, char* argv[]) {
    double h1 = 1e-3;
    double h2 = h1/2;
    double I1 = integr(h1);
    double I2 = integr(h2);
    double E = abs(I1-I2)/(1-0.5);
    for (int i=0; i<30; i++) {
        h1 /= 2;
        h2 = h1/2;
        I1 = integr(h1);
        I2 = integr(h2);
        E = abs(I1-I2)/(1-0.5);
        cout << "Ошибка по правилу Рунге = " << E << endl;
    }

    //определяю прмиерный максимум производной, чтобы оценить нужный шаг по сетке
    /*double H = 3./100;
    double max = -1e5;
    for (int i=0; i<100; i++){
        double x = i*H;
        if (dev2(x) > max) {
            max = dev2(x);
        }
    }*/
    // max|f''| = 7500

    return 0;
}
