#include <iostream>
#include <vector>
#include <cmath>
#include "run_through_3d.h"
#include <fstream>

using namespace std;


class spline_3 {
public:
    spline_3 () {}
    spline_3 (const vector<double>& x_pattern, const vector<double>& y_pattern) {
        nodes = x_pattern;
        value = y_pattern;

        n = nodes.size();
        vector<double> v1 = {0};
        vector<double> v2 (n-2, 2);
        vector<double> v3 = {(nodes[2]-nodes[1])/(nodes[2]-nodes[0])};
        vector<double> v4 = {6*(( (value[2]-value[1])/(nodes[2]-nodes[1]) - (value[1]-value[0])/(nodes[1]-nodes[0]) ) / (nodes[2]-nodes[0]))};
        for (int i=2; i<n-2; i++) {
            v1.push_back((nodes[i]-nodes[i-1])/(nodes[i+1]-nodes[i-1]));
            v3.push_back((nodes[i+1]-nodes[i])/(nodes[i+1]-nodes[i-1]));
            v4.push_back(6*( ((value[i+1]-value[i])/(nodes[i+1]-nodes[i]) - (value[i]-value[i-1])/(nodes[i]-nodes[i-1])) / (nodes[i+1]-nodes[i-1])) );
        }
        v1.push_back((nodes[n-2]-nodes[n-3])/(nodes[n-1]-nodes[n-3]));
        v3.push_back(0);
        v4.push_back(6*(( (value[n-1]-value[n-2])/(nodes[n-1]-nodes[n-2]) - (value[n-2]-value[n-3])/(nodes[n-2]-nodes[n-3]) ) / (nodes[n-1]-nodes[n-3])));
        run_through coef_c (v1, v2, v3, v4);
        c = coef_c.get_solution();
        // добавим 0 в конец вектора с, тк у нас вторая производная последнего сплйна равна 0 (мы так задали)
        c.push_back(0);
        b.push_back(c[0]*(nodes[1]-nodes[0])/3 + (value[1]-value[0])/(nodes[1]-nodes[0]));
        d.push_back(c[0]/(nodes[1]-nodes[0]));
        a.push_back(value[1]);
        for (int i=1; i<n-1; i++) {
            a.push_back(value[i+1]);
            b.push_back(c[i]*(nodes[i+1]-nodes[i])/3 + c[i-1]*(nodes[i+1]-nodes[i])/6 + (value[i+1]-value[i])/(nodes[i+1]-nodes[i]));
            d.push_back((c[i]-c[i-1])/(nodes[i+1]-nodes[i]));
        }
    }
    double get_value (double x) {
        double s=0;
        for (int i=1; i<n; i++) {
            if (x>=nodes[i-1] && x<=nodes[i]) {
                s = a[i-1] + b[i-1]*(x-nodes[i]) + c[i-1]/2*pow(x-nodes[i], 2)
                  + d[i-1]/6*pow(x-nodes[i], 3);
                break;
            }
        }
        return s;
    }
    double derivative (double x) {
        double s=0;
        for (int i=1; i<n; i++) {
            if (x>=nodes[i-1] && x<=nodes[i]) {
                s = b[i-1] + c[i-1]*pow(x-nodes[i], 1)
                  + d[i-1]/2*pow(x-nodes[i], 2);
                break;
            }
        }
        return s;
    }

private:
    vector<double> a , b, c, d;
    vector<double> nodes;
    vector<double> value;
    int n;
};

double der2 () {
	vector<double> x_base1 = {0.1, 0.3, 0.42, 0.5, 0.73};
	vector<double> y_base1 = {0, -3.1, 1.6, 1.45, 2.34};
	spline_3 N1(x_base1, y_base1);
	double h = 1e-5;
	double max=0;
	double H = 0.63/1000;
	for (int i=0; i<1001; i++) {
		if (abs((N1.get_value(H*i+h)-2*N1.get_value(H*i)+N1.get_value(H*i-h))/(h*h)) >= max) {
			max = abs((N1.get_value(H*i-h)-2*N1.get_value(H*i)+N1.get_value(H*i+h))/(h*h));
		}
	}
	return max;
}

int main(int argc, char* argv[]) {
	vector<double> x_base = {0.1, 0.3, 0.42, 0.5, 0.73};
	vector<double> y_base = {0, -3.1, 1.6, 1.45, 2.34};

	spline_3 s(x_base, y_base);
    cout << s.derivative(0.47) << endl;
    cout << (s.get_value(0.47+1e-5)-s.get_value(0.47))/(1e-5);
    cout << endl << der2();
    return 0;
}
