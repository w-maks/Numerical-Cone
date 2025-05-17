#include "rk4.h"
#include "cone.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>


int main() {
    const std::vector<std::vector<double>> wp = {{1.1, 1.0, 0., 0.}, {3., 0.5, 2., 1.}}; //wp[i] {phi0, z0, omega0, vz0}
    constexpr size_t n = 4;
    std::vector<double> s(n), a(n);
    std::array<double, 3> xyz{};

    std::function<void(double, const std::vector<double>&, std::vector<double>&)> f = cone_derivatives;

    const std::vector<std::vector<std::string>> files = {{"cone1.csv","xyz1.csv"}, {"cone2.csv", "xyz2.csv"}};

    for(int i = 0; i < 2; i++) {
        constexpr int N = 500;
        double t = 0.0;
        s[0] = wp[i][0];
        s[1] = wp[i][1];
        s[2] = wp[i][2];
        s[3] = wp[i][3];

        std::ofstream fc(files[i][0]);
        fc << std::fixed << std::setprecision(6);
        fc << "t,phi,z,omega,vz,E,E_anal\n";

        std::ofstream fx(files[i][1]);
        fx << std::fixed << std::setprecision(6);
        fx << "t,x,y,z\n";

        for (int j = 0; j <= N; ++j) {
            constexpr double dt = 0.1;
            fc << t << "," << s[0] << "," << s[1] << "," << s[2] << "," << s[3] << "," << E(T(s[1], s[2], s[3]), U(s[0], s[1])) << "," << E_anal(wp[i]) << "\n";
            xyz = transform_cone_to_lab(s[0], s[1]);
            fx << t << "," << xyz[0] << "," << xyz[1] << "," << xyz[2] << "\n";

            rk4_vec(t, dt, s, f);
            t += dt;
        }
    }
    return 0;
}
