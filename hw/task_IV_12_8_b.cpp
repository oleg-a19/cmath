#include <iostream>
#include <cmath>

using namespace std;


int main(int argc, char* argv[]) {
    // Ищем первый корень
    double x0 = 1/sqrt(2);
    double x = 0.5/sqrt(2)*exp(x0*x0-0.5);
    while (abs(x-x0) > 0.5*1e-3) {
        x0 = x;
        x = 0.5/sqrt(2)*exp(x0*x0-0.5);
    }

    // Ищем второй корень
    double x02 = 1/sqrt(2);
    double x2 = sqrt(log(2*sqrt(2)*x02)+0.5);
    while (abs(x2-x02) > 0.5*1e-3) {
        x02 = x2;
        x2 = sqrt(log(2*sqrt(2)*x02)+0.5);
    }
    cout << fixed;
    cout.precision(6);
    cout << "x1 = " << x << "\nx2 = " << x2
        << "\nШирина функции на полувысоте l = " << x2-x;
    return 0;
}
