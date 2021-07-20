#include <string>
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include <unordered_map>

#include "bed_utils.hpp"

using namespace std ;

// list of all loaded tracks
vector<Track> bed_tracks ;

string get_svtype(uint8_t svtype) {
    if (svtype == SVTYPE_DEL) {
        return "DEL" ;
    }
    if (svtype == SVTYPE_INS) {
        return "INS" ;
    }
    return "INV" ;
}

uint8_t parse_svtype(string svtype) {
    if (svtype == "DEL") {
        return SVTYPE_DEL ;
    }
    if (svtype == "INS") {
        return SVTYPE_INS ;
    }
    return SVTYPE_MISC ;
}


uint8_t get_chromosome_index(string chrom) {
    string ch = chrom.substr(3, chrom.length() - 3) ;
    if (ch == "X") {
        return 22 ;
    }
    if (ch == "Y") {
        return 23 ;
    }
    return std::stoi(ch) ;
}

string get_chromosome_name(uint8_t index) {
    if (index < 22) {
        return "chr" + std::to_string(index) ;
    } else if (index == 22) {
        return "chrX" ;
    } else if (index == 23) {
        return "chrY" ;
    } else {
        cout << int(index) << endl ;
        assert(false) ;
        return "chrUn" ;
    }
}

string SimpleTrack::get_name() const {
    std::string base = get_chromosome_name(chrom) + "_" + std::to_string(begin) + "_" + std::to_string(end) ;
    if (svtype == SVTYPE_DEL) {
        return "DEL@" + base ;
    }
    if (svtype == SVTYPE_INS) {
        return "INS@" + base ;
    }
    return base ;
}

void SimpleTrack::parse_from_name(string name) {
    stringstream ss(name) ;
    vector<string> tokens ;
    string token ;
    while (getline(ss, token, '_')) {
        tokens.push_back(token) ;
    }
    int i = tokens[0].find('@') ;
    svtype = parse_svtype(tokens[0].substr(0, 3)) ;
    chrom = get_chromosome_index(tokens[0].substr(i + 1, tokens[0].length() - (i + 1))) ;
    begin = stoi(tokens[1]) ;
    end = stoi(tokens[2]) ;
}

std::vector<Track> load_tracks_from_file(string path) {
    cout << "Parsing BED file: " << path << ".." << endl ;
    std::ifstream bed_file(path) ;
    std::string line ;
    int i = 0 ;
    std::vector<Track> tracks ;
    unordered_map<string, int> header ;
    while (std::getline(bed_file, line)) {
        istringstream iss(line) ;
        if (i == 0) {
            if (line[0] != '#') {
                cout << "BAM header not present (first line doesn't begin with #). Aborting.." << endl ;
                exit(0) ;
            }
            vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}} ;
            for (int i = 0; i < tokens.size(); i++) {
                header[tokens[i]] = i ;
            }
            i += 1 ;
        } else {
            vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}} ;
            Track track ;
            // make sure all chromosome names begin with "chr"
            string _chrom = tokens[0] ;
            if (_chrom.substr(0, 3) != "chr") {
                _chrom = "chr" + _chrom ;
            } 
            track.chrom = get_chromosome_index(_chrom) ;
            track.begin = std::stoi(tokens[1]) ;
            track.end = std::stoi(tokens[2]) ;
            if (header.find("SEQ") != header.end()) {
                track.seq = tokens[header["SEQ"]] ;
            }
            if (header.find("GENOTYPE") != header.end()) {
                track.genotype = tokens[header["GENOTYPE"]] ;
            }
            if (header.find("SVTYPE") != header.end()) {
                track.svtype = parse_svtype(tokens[header["SVTYPE"]]) ;
            } else {
                track.svtype = SVTYPE_DEL ; 
            }
            if (header.find("SVLEN") != header.end()) {
                track.svlen = abs(std::stoi(tokens[header["SVLEN"]])) ;
            } else {
                // calculate SVLEN based on other parameters
                if (track.svtype != SVTYPE_INS) {
                    track.svlen = track.end - track.begin ;
                } else {
                    track.svlen = track.seq.length() ;
                }
            }
            tracks.push_back(track) ;
        }
    }
    cout << "Loaded " << tracks.size() << " tracks." << endl ;
    return tracks ;
}

std::unordered_map<Track, int> load_tracks_from_file_as_dict(string path) {
    auto _tracks = load_tracks_from_file(path) ;
    std::unordered_map<Track, int> tracks ;
    for (auto track = _tracks.begin(); track != _tracks.end(); track++) {
        tracks[*track] = 1 ;
    }
    return tracks ;
}

