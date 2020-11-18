#include <string>
#include <iostream>

#include "config.hpp"

using namespace std ;

Configuration* Configuration::instance = nullptr ;

Configuration* Configuration::getInstance() {
    if (instance == nullptr) {
        instance = new Configuration() ;
    }
    return instance ; 
}

Configuration::Configuration() :
    parser("Nebula", "Ultra-efficient mapping-free structural variation genotyper.") {
    parser.add_options()
        ("cgc", "", cxxopts::value<bool>())
        ("bed", "", cxxopts::value<std::vector<std::string>>())
        ("bam", "", cxxopts::value<std::vector<std::string>>())
        ("std", "", cxxopts::value<int>())
        ("fastq", "", cxxopts::value<std::vector<std::string>>())
        ("reduce", "", cxxopts::value<bool>())
        ("unique", "", cxxopts::value<bool>())
        ("select", "", cxxopts::value<bool>())
        ("threads", "", cxxopts::value<int>())
        ("workdir", "", cxxopts::value<std::string>())
        ("kmers", "", cxxopts::value<std::vector<std::string>>())
        ("samples", "", cxxopts::value<std::vector<std::string>>())
        ("gc_kmers", "", cxxopts::value<std::string>())
        ("reference", "", cxxopts::value<std::string>())
        ("depth_kmers", "", cxxopts::value<std::string>())
    ;
}

void Configuration::parse(int argc, char** argv) {
    auto results = parser.parse(argc, argv) ;

    cgc = false ;
    if (results.count("cgc")) {
        cgc = true ;
    }
    if (results.count("bed")) {
        bed = results["bed"].as<std::vector<std::string>>() ;
    }
    if (results.count("bam")) {
        bam = results["bam"].as<std::vector<std::string>>() ;
    }
    if (results.count("fastq")) {
        fastq = results["fastq"].as<std::vector<std::string>>() ;
    }
    if (results.count("reduce")) {
        reduce = true ; 
    } else {
        reduce = false ;
    }
    if (results.count("unique")) {
        unique = true ; 
    } else {
        unique = false ;
    }
    if (results.count("select")) {
        select = true ; 
    } else {
        select = false ;
    }
    if (results.count("threads")) {
        threads = results["threads"].as<int>() ;
    }
    if (results.count("workdir")) {
        workdir = results["workdir"].as<std::string>() ;
    }
    if (results.count("kmers")) {
        kmers = results["kmers"].as<std::vector<std::string>>() ;
    }
    if (results.count("samples")) {
        samples = results["samples"].as<std::vector<std::string>>() ;
    }
    if (results.count("gc_kmers")) {
        gc_kmers = results["gc_kmers"].as<string>() ;
    }
    if (results.count("depth_kmers")) {
        depth_kmers = results["depth_kmers"].as<string>() ;
    }
    if (results.count("reference")) {
        reference = results["reference"].as<string>() ;
    }
}

