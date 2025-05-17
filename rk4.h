#ifndef RK4_H
#define RK4_H

#include <vector>
#include <functional>

void rk4_vec(double t, double dt, std::vector<double>& s, const std::function<void(double, const std::vector<double>&, std::vector<double>&)>& f);

#endif
