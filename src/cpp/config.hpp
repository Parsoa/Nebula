#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <string>
#include <vector>

#include "cxxopts.hpp"
#include "bed_utils.hpp"

class Configuration {

private:
    static Configuration* instance ;

public:
    static Configuration* getInstance() ;

    void parse(int argc, char* argv[]) ;

    bool cgc ;
    bool reduce ;
    bool unique ;
    bool select ;
    int threads ;
    std::string workdir ;
    std::string gc_kmers ;
    std::string depth_kmers ;
    std::string reference ;
    std::vector<std::string> bed ;
    std::vector<std::string> bam ;
    std::vector<std::string> fastq ;
    std::vector<std::string> kmers ;
    std::vector<std::string> samples ;
    // map of all loaded tracks
    std::unordered_map<SimpleTrack, Track> tracks ;

private:

    Configuration() ;

    Configuration(Configuration const&) = delete ;
    void operator=(Configuration const&) = delete ;

    Configuration& operator[](std::string) ;
    
    cxxopts::Options parser ;
};

#endif
