#include <iomanip>

#include "chromosomes.hpp"

using namespace std ;

std::vector<std::string> chromosomes ;
std::unordered_map<std::string, char*> chromosome_seqs ;

int get_reference_size(ifstream &fasta_file) {
    fasta_file.seekg(0, ios_base::end) ;
    int l = fasta_file.tellg() ;
    fasta_file.seekg(0, ios_base::beg) ;
    return l ;
}

void load_chromosomes(string path) {
    cout << "Loading reference genome.." << endl ;
    ifstream fasta_file ;
    fasta_file.open(path, ios::binary) ;
    // maximum size of a chromosome, kinda arbitrary
    char* buffer = (char*) malloc(sizeof(char) * 300000000) ;
    int state ;
    uint64_t n = 0 ;
    std::string line ;
    std::getline(fasta_file, line) ;
    while (true) {
        if (line[0] == '>') {
            int o = 3 ; //hg19 style
            if (line.substr(1, 3) == "chr") {
                o = 0 ;
            }
            int l = line.length() ;
            if (l == 6 - o || l == 5 - o) {
                if (line[4 - o] == 'X' || line[4 - o] == 'Y' || (line[4 - o] >= '1' && line[4 - o] <= '9')) {
                    string chrom = (o == 3 ? "chr" : "") + line.substr(1, l - 1) ; // chromosome names will begin with "chr"
                    //cout << "Collecting " << chrom << ".." << endl ;
                    while(std::getline(fasta_file, line)) {
                         if (line[0] == '>') {
                             break ;
                         }
                         for (int i = 0; i < line.length(); i++) {
                            line[i] = toupper(line[i]) ;
                         }
                         memcpy(buffer + n, line.c_str(), line.length()) ;
                         n += line.length() ;
                    }
                    buffer[n] = '\0' ;
                    cout << "Extracted " << std::setw(6) << std::left << chrom << " with " << std::setw(11) << left << n << " bases." << endl ;
                    char* s = (char*) malloc(sizeof(char) * (n + 1)) ;
                    memcpy(s, buffer, n + 1) ;
                    chromosomes.push_back(chrom) ;
                    chromosome_seqs[chrom] = s ;
                    n = 0 ;
                    continue ;
                }
            }
        } 
        if (!std::getline(fasta_file, line)) {
            break ;
        }
    }
    free(buffer) ;
}
