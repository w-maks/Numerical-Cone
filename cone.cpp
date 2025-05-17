#include "cone.h"

using namespace constants;
void cone_derivatives(double t, const std::vector<double>& s, std::vector<double>& k)
{
    k[0] = s[2];
    k[1] = s[3];
    k[2] = -g * std::pow(std::cos(alfa), 2) / std::sin(alfa) * std::sin(s[0]) / s[1] - 2.0 * s[2] * s[3] / s[1];
    k[3] = std::pow(std::sin(alfa), 2) * s[1] * std::pow(s[2], 2) - g * std::sin(alfa) * std::pow(std::cos(alfa), 2) * (1 - std::cos(s[0]));
}

double T(const double z, const double omega, const double vz) {
    return 0.5 * ( std::pow(z * std::tan(alfa) * omega, 2) + std::pow(vz / std::cos(alfa), 2));
}

double U(const double phi, const double z) {
    return g * z * std::sin(alfa) * (1 - std::cos(phi));
}

double E(const double T, const double U) {
    return T + U;
}

double E_anal(const std::vector<double> &wp) {
    return T(wp[1],wp[2],wp[3]) + U(wp[0],wp[1]) ;
}

std::array<double, 3> transform_cone_to_lab(const double phi, const double z) {
    const double rho = z * std::tan(alfa);
    const double x = rho * std::cos(phi);
    const double y = rho * std::sin(phi);

    constexpr double theta = M_PI / 2.0 - alfa;
    const double sin_theta = std::sin(theta);
    const double cos_theta = std::cos(theta);

    const double x_prime =  cos_theta * x + sin_theta * z;
    const double y_prime =  y;
    const double z_prime = -sin_theta * x + cos_theta * z;

    return {x_prime, y_prime, z_prime};
}
