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
        vector<double> v2 (n-2, 2); // вот тут лажа я думаю, так как размер должен быть n, как у остальных
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

private:
    vector<double> a , b, c, d;
    vector<double> nodes;
    vector<double> value;
    int n;
};

int main(int argc, char* argv[]) {
    /*cout << "проверка прогонки";
    vector<double> a1,b1,c1,d1;
    a1={0,2,-1,8}; b1={1,-2,-4,1}; c1={-6,4,6,0}; d1={45,-36,3,-79};
    run_through el(a1,b1,c1,d1);
    vector<double> x = el.get_solution();
    for (const auto& i : x) cout << i << " ";
    cout << "\nend of check";*/

    double a = 0.;
	double b = 2*M_PI;
	int n = 5;
	vector<double> x_base(n, a);

	for (int i=1; i<n; i++) {
		x_base[i] += i*(b-a)/(n-1);
	}
	vector<double> y_base(n);
	for (unsigned int i=0; i<x_base.size(); i++) {
        y_base[i] = sin(x_base[i]);
	}

	spline_3 s(x_base, y_base);
    cout << s.get_value(M_PI/6);
    /*vector<double> test_x(20, a+0.1);
    for (int i=0; i<20; i++) {
        test_x[i]+= i*(b-(a+0.1))/(20-1);
    }
	vector<double> test_y(test_x.size());
    ofstream outf("data1.txt");
	for (unsigned int i=0; i<test_x.size(); i++) {
		test_y[i] = s.get_value(test_x[i]);
        outf << test_x[i] << " " << test_y[i] << endl;
	}*/
    return 0;
}
