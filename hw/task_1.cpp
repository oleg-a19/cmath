#include <iostream>
#include <cmath>

using namespace std;


int main(int argc, char* argv[]) {
    // Ищем первый корень
    double y0 = 0.5;
    double x0 = 0.5;
    double x = atan(y0);
    double y = sqrt(1-x0*x0);
    while (abs(y-y0) > 1e-6 && abs(x-x0) > 1e-6) {
        y0 = y;
        x0 = x;
        x = atan(y0);
        y = sqrt(1-x0*x0);
    }

    // Ищем второй корень
    double y02 = -0.5;
    double x02 = -0.5;
    double x2 = atan(y02);
    double y2 = -sqrt(1-x02*x02);
    while (abs(y2-y02) > 1e-6 && abs(x2-x02) > 1e-6) {
        y02 = y2;
        x02 = x2;
        x2 = atan(y02);
        y2 = -sqrt(1-x02*x02);
    }
    cout << fixed;
    cout.precision(6);
    cout << "x1 = " << x << ", y1 = " << y
        << "\nx2 = " << x2 << ", y2 = " << y2;
    return 0;
}
