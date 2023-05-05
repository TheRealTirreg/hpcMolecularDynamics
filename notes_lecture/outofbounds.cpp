#include "Eigen/Dense"

double square_and_sum(Eigen::ArrayXd v, int n) {
    double res = 0;
    for (int i = 1; i < n; ++i) {
        v[i] *= v[i];
        res += v[i];
    }
    return res;
}

int main() {
    int n = 5;
    Eigen::ArrayXd v1(n);
    square_and_sum(v1, n);
    return 0;
}
