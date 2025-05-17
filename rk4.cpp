#include "rk4.h"

void rk4_vec(const double t, const double dt, std::vector<double>& s, const std::function<void(double, const std::vector<double>&, std::vector<double>&)>& f)
{
    const size_t n = s.size();
    std::vector<double> k1(n), k2(n), k3(n), k4(n), w(n);

    w = s;
    f(t, w, k1);

    for (int i = 0; i < n; ++i) w[i] = s[i] + 0.5 * dt * k1[i];
    f(t + 0.5 * dt, w, k2);

    for (int i = 0; i < n; ++i) w[i] = s[i] + 0.5 * dt * k2[i];
    f(t + 0.5 * dt, w, k3);

    for (int i = 0; i < n; ++i) w[i] = s[i] + dt * k3[i];
    f(t + dt, w, k4);

    for (int i = 0; i < n; ++i)  s[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
}
