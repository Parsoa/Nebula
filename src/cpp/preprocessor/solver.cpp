#include "solver.hpp"

void LpSolver::estimate_coverage() {
    std::vector<int> counts ;
    for (auto kmer = kmers.begin(); kmer != kmers.end(); kmer++) {
        if (kmer->second.type == KMER_TYPE_DEPTH) {
            counts.push_back(kmer->second.count) ;
        }
    }
}
