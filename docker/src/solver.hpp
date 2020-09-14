#ifndef SLV_HPP
#define SLV_HPP

#include <vector>
#include <string>
#include <iostream>
#include <iterator>

#include "kmer_utils.hpp"

class LpSolver {

public:


private:

void estimate_coverage() ;

double std ;
double mean ;

std::unordered_map<uint64_t, Kmer> kmers ;

} ;

#endif
