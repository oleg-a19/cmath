#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>

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

double der2 () {
	vector<double> x_base1 = {0.1, 0.3, 0.42, 0.5, 0.73};
	vector<double> y_base1 = {0, -3.1, 1.6, 1.45, 2.34};
	interpolator3_newton N1(x_base1[0], x_base1[x_base1.size()-1], x_base1, y_base1);
	double h = 1e-5;
	double max=0;
	double H = 0.63/1000;
	for (int i=0; i<1001; i++) {
		if (abs((N1.print_y(H*i+h)-2*N1.print_y(H*i)+N1.print_y(H*i-h))/(h*h)) >= max) {
			max = abs((N1.print_y(H*i-h)-2*N1.print_y(H*i)+N1.print_y(H*i+h))/(h*h));
		}
	}
	return max;
}

int main(int argc, char* argv[]) {
	int n=3;
	//vector<double> x_base = {-3, 0, 4};
	//vector<double> y_base = {exp(-3), exp(0), exp(4)};
	vector<double> x_base(n);
	vector<double> y_base(n);
	for (int i=0; i<n; i++) {
		x_base[i] = 0.5*(1+7*cos((M_PI/2+M_PI*i)/n));
		cout << x_base[i] << endl;
		y_base[i] = exp(x_base[i]);
	}
	cout << endl;
	
	/*reverse(x_base.begin(), x_base.end());
	for (int i=0; i<n; i++) {
		cout << x_base[i] << endl;
		y_base[i] = exp(x_base[i]);
	}*/

	interpolator3_newton N(x_base[0], x_base[x_base.size()-1], x_base, y_base);
	double y = N.print_y(1);
	cout << "Answer "<< y;
	N.print_coef();

	/*vector<double> x_base1 = {0.1, 0.3, 0.42, 0.5, 0.73};
	vector<double> y_base1 = {0, -3.1, 1.6, 1.45, 2.34};
	interpolator3_newton N1(x_base1[0], x_base1[x_base.size()-1], x_base1, y_base1);
	double h = 1e-2;
	double m = der2();
	cout << "max 2der " << m << endl;
	vector<double> x_base1 = {0.1, 0.3, 0.42, 0.5, 0.73};
	vector<double> y_base1 = {0, -3.1, 1.6, 1.45, 2.34};
	interpolator3_newton N1(x_base1[0], x_base1[x_base1.size()-1], x_base1, y_base1);
	double h = 1e-5;
	double y1 = (N1.print_y(0.47+h)-N1.print_y(0.47-h))/(2*h);
	cout << y1;*/
	return 0;
}
