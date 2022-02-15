#ifndef RUN_THROUGH_H
#define RUN_THROUGH_H
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;


class run_through {
public:
    run_through () {}
    run_through (const vector<double>& a, const vector<double>& b,
                 const vector<double>& c, const vector<double>& d) {
        int n = a.size();
        x.resize(n);
        vector<double> p,q;
        p.push_back(-c[0]/b[0]);
        q.push_back(d[0]/b[0]);
        // рассчитаем прогоночные коэффициенты
        for (int i=1; i<n-1; i++) {
            p.push_back(-c[i] / ( a[i]*p[i-1]+b[i] ));
            q.push_back(( d[i]-a[i]*q[i-1] ) / ( a[i]*p[i-1]+b[i] ));
        }
        cout << "должен быть от 0 до 1 \np1 = "<<p[1] << endl;
        x[n-1] = ( d[n-1]-a[n-1]*q[n-2] ) / ( a[n-1]*p[n-2]+b[n-1] );
        for (int i=n-2; i>=0; i--) {
            x[i] = p[i]*x[i+1]+q[i];
        }
    }
    vector<double> get_solution () {
        return x;
    }
private:
    vector<double> x;
};

#endif
