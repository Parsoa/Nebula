#include <fstream>
#include <iostream>

#include "config.hpp"
#include "logger.hpp"

using namespace std ;

void print_kmer(Kmer& kmer) {
    nlohmann::json payload ;
    cout << "\"" << decode_kmer(kmer.seq) << "\":" ;
    cout << nlohmann::json(kmer).dump(4) ;
    cout << endl ;
}
