#include <math.h>

#include "stats.hpp"

double NormalDistribution::prob(double x) {
    //TODO
    // a·=·math.exp(-1·*·((x·-·self.mean)·**·2)·/·(2.0·*·self.var))
    // return·(1.0·/·math.sqrt(2.0·*·math.pi·*·self.var))·*·a
    return x ;
}

double NormalDistribution::log_prob(double x) {
    //return·(-1·*·(math.log(self.std)·+·0.5·*·(math.log(2)·+·math.log(math.pi))))·-·(0.5·*·((x·-·self.mean)·**·2·/·(self.std)))
    //
    return (-1 * (log(std) + 0.5 * (log(2) + log(M_PI)))) - (0.5 * ((x - mean) * (x - mean) / std)) ;
}
