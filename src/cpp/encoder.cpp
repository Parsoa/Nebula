#include "encoder.hpp"

#include <cstdint>
#include <math.h>

uint8_t* ReadEncoder::encode(std::string read) {
    int l = read.length() ;
    this->bead = new uint8_t[int(ceil(l / 4))] ;
    int j = 0 ;
    uint8_t* b = this->bead ;
    *b = 0 ;
    for (int i = 0; i < l; i++) {
        char r = read[i] ;
        uint8_t c = 0 ? r == 'A' : 1 ? r == 'C' : 2 ? r == 'G' : 3 ;
        this->bead[i] &= r << ((3 - j) * 1) ;
        j++ ;
        if (j >> 2 == 0) {
            j = 0 ;
            b++ ;
            *b = 0 ;
        }
    }
    return this->bead ;
}

void ReadEncoder::decode() {

}