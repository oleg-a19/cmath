#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;


class interpolator3_newton {
public:
	interpolator3_newton () {}
	interpolator3_newton (const double& begin_segment, const double& end_segment,
				  const vector<double>& x_pattern, const vector<double>& y_pattern) {
		n = x_pattern.size();
		a = begin_segment;
		b = end_segment;
		nodes = x_pattern;
		value = y_pattern;

		// найдём разделённые разности
		vector<double> new_value(n-1);
		vector<double> tmp = value;
		coef.push_back(value[0]);
		for (int i=0; i<n-1; i++) {
			new_value.resize(n-i-1);
			for (int k=0; k<n-i-1; k++) {
				new_value[k] = (tmp[k+1]-tmp[k])/(nodes[i+k+1]-nodes[k]);
			}
			tmp.resize(new_value.size());
			tmp = new_value;
			coef.push_back(new_value[0]); // сохраняем нужную разделённую разность
		}
	}
	// Находим значение интерполяционного многочлена в заданной точке
	double print_y (const double& x) {
		double y = 0;
		for (int i=0; i<n; i++) {
			double P = 1; // произведение скобок (x - x_i) i=0,1..n
			for (int j=0; j<i; j++) {
				P *= (x-nodes[j]);
			}
			y += coef[i]*P;
		}
		return y;
	}
	void print_coef () {
		cout << endl;
		for (int i=0; i<coef.size(); i++) {
			cout << "a"<<i<<" = " << coef[i] << endl;
		}
	}

private:
	int n;
	double a, b;
	vector<double> nodes;
	vector<double> value;
	vector<double> coef; // коэффициенты интерполяционного многочлена
};


int main(int argc, char* argv[]) {
	vector<double> x_base = 

	vector<double> y_base =


	interpolator3_newton N(a, b, x_base, y_base);
	vector<double> test_x = {a+1, (b-a)/2, b-1};
	vector<double> test_y(test_x.size());
	for (unsigned int i=0; i<test_x.size(); i++) {
		test_y[i] = N.print_y(test_x[i]);
		cout << test_x[i] << " " << test_y[i] << endl;
	}
	N.print_coef();

	return 0;
}
