#include <math.h>

#include "stats.hpp"

double NormalDistribution::prob(double x) {
    //double a = math.exp(-1 * ((x - mean) * (x - mean)) / (2.0 * var)) ;
    //return (1.0 / sqrt(2.0 * math.pi * self.var)) * a
    return 0.0 ;
}

double NormalDistribution::log_prob(double x) {
    return (-1 * (log(std) + 0.5 * (log(2) + log(M_PI)))) - (0.5 * ((x - mean) * (x - mean) / std)) ;
}
