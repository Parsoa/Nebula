#include <iostream>
#include <string>
#include "kmc_file.h"

using namespace std ;

CKMCFile kmer_database ;
bool loaded = false ;

void load_database() {
    kmer_database.OpenForRA("/share/hormozdiarilab/Experiments/KMC/HG00512/mer_counts") ;
    loaded = true ;
}

int get_kmer_count(string kmer) {
    if (!loaded){
        load_database() ;
    }
    uint32 count = 0 ;
    CKmerAPI* kmer_object = new CKmerAPI(31) ;
    kmer_object->from_string(kmer) ;
    kmer_database.CheckKmer(*kmer_object, count) ;
    return int(count) ; 
}

int main() {
    cout << get_kmer_count("GCAAACATAGTGAAACCCCGTCTCTACTAAA") ;
    return 0 ;
}
