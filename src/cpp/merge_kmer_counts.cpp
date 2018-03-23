#include <memory>
#include <fstream>
#include <iostream>

#include "json.hpp"

using namespace std ;
using json = nlohmann::json ;

int main() {

    string path = "/share/hormozdiarilab/Codes/NebulousSerendipity/output/CHM1_Lumpy.Del.100bp.DEL.bed/31/CountKmersExactJob" ;

    json kmers ;
    std::ifstream i(path + "/batch_0.json");
    i >> kmers ;
    cout << "total kmers: " << kmers.size() << endl ;

    for (int i = 1 ; i < 48 ; i++) {
        cout << "reading batch " << i << endl ;
        std::auto_ptr<json> batch(new json()) ;
        std::ifstream f(path + "/batch_" + to_string(i) + ".json");
        f >> *batch ;
        for (json::iterator it = batch->begin(); it != batch->end(); ++it) {
            kmers[it.key()] = int(kmers[it.key()]) + int(it.value()) ;
        }
    }

    cout << "done mergin counts, dumping..." << endl ;
    std::ofstream o(path + "/kmers.json") ;
    o << std::setw(4) << kmers << std::endl ;
}
