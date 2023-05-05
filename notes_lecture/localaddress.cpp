#include "Eigen/Dense"

/* double* unit_3d_vector() {
    double v[3] = {0, 0, 1};
    return v;
}*/

Eigen::ArrayXd unit_3d_vector() {
    Eigen::ArrayXd ret(3);
    ret << 0, 0, 1;
    return ret;
}

// solution:
// 1) use vector (but is slower because memory is stored on heap)
// 2) use array (equivalent)
// 3) use Eigen

int main() {
    Eigen::ArrayXd n = unit_3d_vector();
    for (int i = 0; i < 3; ++i)
        n[i] += 1;
    return 0;
}