#include <iostream>
#include <vector>
#include <cmath>

#include <fstream>

std::vector<double> f (double x, double y) {
    double m = 1000.;
    std::vector<double> func(2);
    //func[0] = x-2;
    //func[1] = x+2*y-3;

    func[0] = y;
    func[1] = m*(1-y*y)*y-x;

    return func;
}

struct solution {
    std::vector<double> x;
    std::vector<double> y; // derivative
    std::vector<double> nodes;
};

class reley {
public:
    reley (double start_value1, double start_value_deriveative1, int time1, double step1) {
        start_value = start_value1;
        start_value_deriveative = start_value_deriveative1;
        time = time1;
        step = step1;
    }
    solution solve () {
        s.nodes.resize(time+1);
        s.x.resize(time+1);
        s.y.resize(time+1);
        s.nodes[0]=0;
        s.x[0]=start_value;
        s.y[0]=start_value_deriveative;
        for (int i=0; i<time+1; i++) {
            s.nodes[i] = i;
        }
        double x1, y1, x2, y2;
        x1 = start_value;
        y1 = start_value_deriveative;
        x2 = 0;
        y2 = 0;
        for (int i=0; i<time+1; i++) {
            //std::cout<<s.x[i-1]<<" "<<s.y[i-1]<<std::endl;
            s.x[i] = x1;
            s.y[i] = y1;
            //x2 = s.x[i-1] + step/2*(f(s.x[i-1], s.y[i-1])[0] + f(x1, y1)[0]);
            //y2 = s.y[i-1] + step/2*(f(s.x[i-1], s.y[i-1])[1] + f(x1, y1)[1]);

            //std::cout<<f(x1, y1)[0]<<std::endl;
            while (abs(x2-x1) > 1e-6 && abs(y2-y1) > 1e-6) {
                //std::cout << "+++"<<std::endl;
                //std::cout<<"x[i-1] = "<<s.x[i-1]<<" f(x[i-1], y[i-1])[0] = "<<f(s.x[i-1], s.y[i-1])[0]<< " f(x1, y1)[0] = "<<f(x1, y1)[0]<<std::endl;
                //std::cout<<f(x1, y1)[0]<<std::endl;
                x2 = s.x[i] + step/2*(f(s.x[i], s.y[i])[0] + f(x1, y1)[0]);
                y2 = s.y[i] + step/2*(f(s.x[i], s.y[i])[1] + f(x1, y1)[1]);
                x1 = x2;
                y1 = y2;
                //std::cout << "x1 = "<<x1<<" x2 = "<<x2 <<" y1 = "<<y1<<" y2 = "<<y2<<" f[0] = "<<f(x1, y1)[0]<<" f[1] = "<<f(x1, y1)[1]
                //<<" |x1-x2| = "<<abs(x2-x1)<<" |y1-y2| = "<<abs(y2-y1)<<std::endl;
                //std::cout << "x1-x2 = "<<x1-x2<<" y1-y2 = "<<y1-y2 <<std::endl;
                //std::cout<<"x[i-1] = "<<s.x[i-1]<<" f(x[i-1], y[i-1])[0] = "<<f(s.x[i-1], s.y[i-1])[0]<< " f(x1, y1)[0] = "<<f(x1, y1)[0]<<std::endl;
                //std::cout<<x2<<" "<<y2<<std::endl;

            }
            //std::cout << "t = "<<s.nodes[i]<<" x = "<<s.x[i]<<std::endl;
        }
        return s;
    }
private:
    solution s;
    double start_value, start_value_deriveative, step;
    int time;
};

int main (int argc, char* argv[]) {
    std::ofstream file ("/home/senseye3/Study/c_math/odu/out_reley.csv");

    reley r(0, 0.001, 1000, 1);
    solution s1 = r.solve(); // шаг один потому что в задаче такая сетка задана
    for (int i=0; i<s1.nodes.size(); i++) {
        std::cout << "t = " << s1.nodes[i] << " x = " << s1.x[i] << std::endl;
    }
    file << "t,x" << std::endl;
    for (int index = 0; index < s1.x.size(); ++index){
        file << s1.nodes[index] << "," << s1.x[index] << std::endl;
    }
    file.close();
    return 0;
}
