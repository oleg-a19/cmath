#include <iostream>
#include <vector>
#include <fstream>

using namespace std;


class interpolator_lagrange {
public:
	interpolator_lagrange () {}
	interpolator_lagrange(const double& begin_segment, const double& end_segment, 
				  const vector<double>& x_pattern, const vector<double>& y_pattern) {
		n = x_pattern.size();
		a = begin_segment;
		b = end_segment;
		nodes = x_pattern;
		vector<double> coef; // вектор произведений в знаменателях коэффициентов a_i в сумме интерполятора
		for (int i=0; i<n; i++) {
			double P = 1; // произведение в знаменателе i-го коэффициента интерполятора
			for (int j=0; j<n; j++) {
				if (nodes[i] == nodes[j]) continue;
				P *= (nodes[i] - nodes[j]);
			}
			coef.push_back(P);
		}
		// найдём коэффициент i-го члена в сумме интерполятора
		factor.resize(n);
		for (int i=0; i<n; i++) {
			factor[i] = y_pattern[i]/coef[i];
		}
	}
	// Находим значение интерполяционного многочлена в заданной точке
	double print_y (const double& x) {
		double y = 0;
		for (int i=0; i<n; i++) {
			double P = 1; // произведение скобок (x - x_i) i=0,1..n
			for (int j=0; j<n; j++) {
				if (i == j) continue;
				P *= (x-nodes[j]);
			}
			y += factor[i]*P;
		}
		return y;
	}
	
private:
	int n;
	double a, b;
	vector<double> nodes;
	vector<double> factor;
};


int main(int argc, char* argv[]) {
	double a = 0.;
	double b = 10.;
	int n = 20;
	vector<double> x_base(n, a);
	for (int i=1; i<n; i++) {
		x_base[i] += i*(b-a)/(n-1);
	}
	vector<double> y_base(n);
	for (unsigned int i=0; i<x_base.size(); i++) {
		y_base[i] = x_base[i]*x_base[i];
	}
		
	interpolator_lagrange l(a, b, x_base, y_base);
	vector<double> test_x = {a+1, (b-a)/2, b-1};
	vector<double> test_y(test_x.size());
	for (unsigned int i=0; i<test_x.size(); i++) {
		test_y[i] = l.print_y(test_x[i]);
		cout << test_x[i] << " " << test_y[i] << endl;
	}
	
	/*ofstream outf("data.txt");
	for (unsigned int i=0; i<x.size(); i++) {
		outf << x[i] << " " << y[i] << endl;
	}*/
	return 0;
}
