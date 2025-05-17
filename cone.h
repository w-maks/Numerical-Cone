#ifndef CONE_H
#define CONE_H

#include <vector>
#include <array>
#include <cmath>

namespace constants {
    inline constexpr double g = 9.81;
    inline constexpr double alfa = 0.5;
}

void cone_derivatives(double t, const std::vector<double>& s, std::vector<double>& k);

double T(double z, double omega, double vz);

double U(double phi, double z);

double E(double T, double U);

double E_anal(const std::vector<double> &wp);

std::array<double, 3> transform_cone_to_lab(double phi, double z);

#endif
