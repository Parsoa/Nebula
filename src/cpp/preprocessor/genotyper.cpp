#include <omp.h>
#include <mutex>
#include <math.h>
#include <string>
#include <stdio.h>
#include <fstream>
#include <iomanip>
#include <assert.h>
#include <iostream>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unordered_map>
#include <bits/stdc++.h>

#include "stats.hpp"
#include "logger.hpp"
#include "counter.hpp"
#include "genotyper.hpp"

#include "json.hpp"
#include "ortools/linear_solver/linear_solver.h"

using namespace std ;
namespace lp = operations_research ;

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

double calc_mean(vector<int> counts) {
    return accumulate(counts.begin(), counts.end(), 0) / counts.size() ;
}

double calc_std(vector<int> counts, double mean) {
    double sum = 0 ;
    for (auto count: counts) {
        sum += (count - mean) * (count - mean) ;
    }
    double var = sum / counts.size() ;
    return sqrt(var) ;
}

std::string genotypes[] = {"0/0", "0/1", "1/1"} ;
double calc_kmer_genotype_likelihood(Kmer& kmer, double count, string genotype, double mean, double std) {
    if (genotype == "1/0") {
        genotype = "0/1" ;
    }
    double c = genotype == "0/0" ? 0.0 : genotype == "0/1" ? 0.5 : 1.0 ;
    double s = genotype == "0/0" ? 0.25 : genotype == "0/1" ? 0.5 : 1.0 ;
    if (kmer.type == KMER_TYPE_JUNCTION || kmer.trend) {
       return NormalDistribution(c * mean, s * std).log_prob(count) ; 
    } else {
        return NormalDistribution((1 - c) * mean, (0.25 / s) * std).log_prob(count) ; 
    }
}

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

void Genotyper::run() {
    auto c = Configuration::getInstance() ;
    load_kmers() ;
    if (c->reduce) {
        cout << "Loading counted kmers.." << endl ;
        load_counts() ;
    } else {
        count_kmers() ;
    }
    estimate_coverage() ;
    adjust_gc_coverage() ;
    if (c->select) {
        filter_kmers() ;
    }
    cluster_kmers() ;
    solve_lp() ;
    export_genotypes() ;
}

void Genotyper::load_kmers() {
    auto c = Configuration::getInstance() ;
    cout << "loading gc kmers.." << endl ;
    ifstream json_file(c->gc_kmers) ;
    nlohmann::json kmers_json ;
    json_file >> kmers_json ;
    for (nlohmann::json::iterator kmer = kmers_json.begin(); kmer != kmers_json.end(); kmer++) {
        uint64_t k = encode_kmer(canonicalize(kmer.key()).c_str()) ;
        Kmer _kmer ;
        _kmer.gc = kmer.value().at("gc").get<int>() ;
        _kmer.seq = k ;
        gc_kmers.emplace(std::make_pair(k, _kmer)) ;
    }
    json_file.close() ;
    cout << "loading depth kmers.." << endl ;
    json_file.open(c->depth_kmers) ;
    nlohmann::json depth_kmers_json ;
    json_file >> depth_kmers_json ;
    for (nlohmann::json::iterator kmer = depth_kmers_json.begin(); kmer != depth_kmers_json.end(); kmer++) {
        uint64_t k = encode_kmer(canonicalize(kmer.key()).c_str()) ;
        Kmer _kmer ;
        _kmer.seq = k ;
        depth_kmers.emplace(std::make_pair(k, _kmer)) ;
    }
    json_file.close() ;
    cout << "loading genotyping kmers.." << endl ;
    int inner ;
    int junction ;
    int num_batches = 1000 ;
    vector<unordered_map<uint64_t, Kmer>> _kmers(num_batches) ;
    #pragma omp parallel for num_threads(c->threads)
    for (int i = 0; i < num_batches; i++) {
        //cout << "batch" << i << endl ;
        int t = i ;//omp_get_thread_num() ;
        ifstream json_file ;
        string path = c->kmers[0] + "/kmers_batch_" + std::to_string(t) + ".json" ;
        struct stat info ;
        if (stat(path.c_str(), &info) != 0) {
            continue ;
        }
        json_file.open(path) ;
        nlohmann::json kmers_json ;
        json_file >> kmers_json ;
        for (nlohmann::json::iterator kmer = kmers_json.begin(); kmer != kmers_json.end(); ++kmer) {
            uint64_t k = encode_kmer(kmer.key()) ;
            _kmers[t][k] = kmer.value().get<Kmer>() ;
            _kmers[t][k].seq = encode_kmer(kmer.key()) ;
            _kmers[t][k].count = 0 ; 
            _kmers[t][k].total = 0 ;
            if (_kmers[t][k].type == KMER_TYPE_JUNCTION) {
                junction += 1 ;
            } else {
                inner += 1 ;
            }
        }
        json_file.close() ;
    }
    cout << "Loaded all kmers batches, merging.." << endl ;
    for (int i = 0; i < num_batches; i++) {
        kmers.insert(_kmers[i].begin(), _kmers[i].end()) ;
        _kmers[i].clear() ;
    }
    cout << "Loaded all kmers:" << endl ;
    cout << "\t" << gc_kmers.size() << " GC kmers." << endl ;
    cout << "\t" << depth_kmers.size() << " depth kmers." << endl ;
    cout << "\t" << kmers.size() << " genotyping kmers." << endl ;
    cout << "\t\t" << inner << " inner kmers." << endl ;
    cout << "\t\t" << junction << " junction kmers." << endl ;
}

void Genotyper::count_kmers() {
    auto counter = KmerCounter(16, &depth_kmers, &gc_kmers, &kmers) ;
    counter.run() ;
}

void Genotyper::load_counts() {
    auto counter = KmerCounter(16, &depth_kmers, &gc_kmers, &kmers) ;
    counter.load_counts() ;
}

void Genotyper::estimate_coverage() {
    std::vector<int> counts ;
    for (auto kmer = depth_kmers.begin(); kmer != depth_kmers.end(); kmer++) {
        counts.push_back(kmer->second.count) ;
    }
    while (true) {
        mean = calc_mean(counts) ;
        std = calc_std(counts, mean) ;
        counts.erase(std::remove_if(counts.begin(), counts.end(),
            [this](const int & o) { return o > this->mean + 3 * this->std ; }), counts.end()) ;
        double m = calc_mean(counts) ;
        cout << "Coverage " << mean << "x with variance " << std << "." << endl ;
        if (abs(mean - m) < 0.5) {
            break ;
        }
    }
    coverage = mean ;
    cout << "Coverage " << mean << "x with variance " << std << "." << endl ;
}

void Genotyper::adjust_gc_coverage() {
    vector<int> counts ;
    vector<vector<int>> gc(101) ;
    for (auto kmer = gc_kmers.begin(); kmer != gc_kmers.end(); kmer++) {
        gc[kmer->second.gc].push_back(kmer->second.count) ;
    }
    for (int i = 0; i <= 100; i++) {
        if (gc[i].size() > 0) {
            while (true) {
                double _mean = calc_mean(gc[i]) ;
                double _std = calc_std(gc[i], _mean) ;
                gc[i].erase(std::remove_if(gc[i].begin(), gc[i].end(),
                    [_mean, _std](const int & o) { return o > _mean + 3 * _std ; }), gc[i].end()) ;
                double m = calc_mean(gc[i]) ;
                if (abs(_mean - m) < 0.5) {
                    gc_coverage.push_back(max(0.2 * coverage, _mean)) ;
                    break ;
                }
            }
        } else {
            gc_coverage.push_back(coverage) ;
        }
        //cout << gc_coverage[gc_coverage.size() - 1] << endl ;
    }
}

void Genotyper::filter_kmers() {
    auto c = Configuration::getInstance() ;
    cout << "Filtering " << kmers.size() << " kmers.." << endl ;
    int e = 0 ;
    int threads = min(c->threads, int(kmers.bucket_count())) ;
    vector<vector<uint64_t>> filters(threads) ;
    cout << "Filtering " << kmers.size() << " kmers using " << threads << " threads.." << endl ;
    std::ofstream j ;
    j.open(c->workdir + "/error.json") ;
    j << "{\n" ;
    //#pragma omp parallel for num_threads(threads)
    for(size_t b = 0; b < kmers.bucket_count(); b++) {
        for(auto it = kmers.begin(b); it != kmers.end(b); it++) {
            auto kmer = it->second ;
            int t = omp_get_thread_num() ;
            int count = kmer.inverse ? kmer.total - kmer.count : kmer.count ;
            int residue = kmer.inverse ? 0 : kmer.loci.size() - kmer.tracks.size() ;
            int coverage = kmer.type == KMER_TYPE_JUNCTION ? this->coverage : gc_coverage[kmer.gc] ;
            int lp_count = min(count, coverage * int(kmer.loci.size())) ;
            if (kmer.loci.size() > 1) {
                filters[t].push_back(it->first) ;
                continue ;
            }
            if (kmer.loci.size() != 1) {
                //cout << decode_kmer(it->first) << endl ;
                //cout << kmer.loci.size() << endl ;
                e += 1 ;
                //if (e != 0) {
                //    j << ",\n" ;
                //}
                //j << "\"" << decode_kmer(it->first) << "\":" ;
                //j << nlohmann::json(kmer).dump(4) ;
            }
            //assert(kmer.loci.size() == 1) ;
            if (kmer.tracks.size() > 1) {
                filters[t].push_back(it->first) ;
                continue ;
            }
            //if (kmer.tracks.size() != 1) {
            //    print_kmer(kmer) ;
            //}
            assert(kmer.tracks.size() == 1) ;
            unordered_map<string, double> likelihoods ;
            double m = -100000000 ;
            string choice ;
            //cout << decode_kmer(it->first) << " " << lp_count << endl ;
            for (auto g: genotypes) {
                likelihoods[g] = calc_kmer_genotype_likelihood(kmer, lp_count, g, coverage, std) ;
                if (likelihoods[g] > m) {
                    m = likelihoods[g] ;
                    choice = g ;
                }
            //cout << g << ": " << likelihoods[g] << endl ; 
            }
            //cout << "Choice: " << choice << ", Correct: " ;
            auto track = kmer.tracks.begin()->first ;
            //if (c->tracks.find(track) != c->tracks.end()) {
            //    cout << c->tracks.find(track)->second.genotype << endl ;
            //} else {
            //    cout << "0/0" << endl ;
            //}
            bool pass = false ;
            if (choice.find("1") != string::npos) {
                if (c->tracks.find(track) != c->tracks.end() && c->tracks.find(track)->second.genotype.find("1") != string::npos) {
                    pass = true ;
                }
            } else {
                if (c->tracks.find(track) == c->tracks.end() || c->tracks.find(track)->second.genotype.find("1") == string::npos) {
                    pass = true ;
                }
            }
            if (!pass) {
                filters[t].push_back(it->first) ;
                continue ;
            }
        }
    }
    j << "\n}\n" ;
    j.close() ;
    for (auto it = filters.begin(); it != filters.end(); it++) {
        for (auto k: *it) {
            kmers.erase(kmers.find(k)) ;
        }
    }
    int inner = 0 ;
    int junction = 0 ;
    for (auto kmer = kmers.begin(); kmer != kmers.end(); kmer++) {
        if (kmer->second.type == KMER_TYPE_INNER) {
            inner += 1 ;
        } else {
            junction += 1 ;
        }
    }
    cout << kmers.size() << " kmers remain after filtering:" << endl ;
    cout << "\tInner kmers: " << inner << endl ;
    cout << "\tJunction kmers: " << junction << endl ;
    //cout << e << " inconsistent kmers." << endl ;
}

void Genotyper::cluster_kmers() {
    int n = 0 ;
    int t = 0 ;
    for (auto kmer = kmers.begin(); kmer != kmers.end(); kmer++) {
        n += 1 ;
        for (auto track = kmer->second.tracks.begin(); track != kmer->second.tracks.end(); track++) {
            tracks[track->first].push_back(kmer->first) ;
            t += 1 ;
        }
    }
    for (auto track = tracks.begin(); track != tracks.end(); track++) {
        cluster_index[track->first] = -1 ;
    }
    cout << "Clustering " << tracks.size() << " tracks with " << n << " kmers.." << endl ;
    for (auto it = tracks.begin(); it != tracks.end(); it++) {
        auto track = it->first ;
        int cluster = cluster_index[track] ;
        for (auto k = tracks[track].begin(); k != tracks[track].end(); k++) {
            for (auto _it = kmers[*k].tracks.begin(); _it != kmers[*k].tracks.end(); _it++) {
                auto _track = _it->first ;
                int _cluster = cluster_index[_track] ;
                if (track != _track) { // look at other tracks only
                    if (_cluster == -1) { // other track not yet clustered
                        if (cluster == -1) { // current track not yet clustered
                            cluster = clusters.size() ;
                            clusters.push_back({track, _track}) ;
                            cluster_index[track] = cluster ;
                            cluster_index[_track] = cluster ;
                        } else {
                            // current track clustered, will add other track later when it is its turn in the loop
                        }
                    } else { // append current track to cluster
                        if (cluster == -1) {
                            cluster = _cluster ;
                            cluster_index[track] = _cluster ;
                            clusters[_cluster].push_back(track) ; 
                        } else {
                            if (cluster != _cluster) {
                                for (auto __track: clusters[_cluster]) {
                                    cluster_index[__track] = cluster ;
                                    clusters[cluster].push_back(__track) ;
                                }
                                clusters[_cluster].clear() ;
                            }
                        }
                    }
                }
            }
        }
        if (cluster == -1) {
            cluster = clusters.size() ;
            clusters.push_back({track}) ;
            cluster_index[track] = cluster ;
        }
    }
    int m = 0 ;
    for (auto it = clusters.begin(); it != clusters.end(); it++) {
        if ((*it).size() > 0) {
            m += 1 ;
        }
    }
    cout << tracks.size() << " tracks with " << clusters.size() << " clusters." << endl ;
    m = 0 ;
    for (auto track = tracks.begin(); track != tracks.end(); track++) {
        m += track->second.size() ;
        genotyping_tracks[track->first].lp_value = 0.0 ;
    }
    cout << m << " kmers in clusters." << endl ;
}

void sample_lp_solve(int batch) {
    auto c = Configuration::getInstance() ;
    lp::MPSolver solver("cluster_" + std::to_string(batch), lp::MPSolver::CLP_LINEAR_PROGRAMMING) ;
    lp::MPObjective* const objective = solver.MutableObjective() ;
    vector<lp::MPVariable*> track_variables ;
    for (int i = 0; i < 100; i++) {
        track_variables.push_back(solver.MakeNumVar(0.0, 1.0, "t_" + std::to_string(i))) ;
    }
    vector<lp::MPVariable*> error_variables ;
    cout << "building lp.." << endl ;
    for (int i = 0; i < 1000; i++) {
        error_variables.push_back(solver.MakeNumVar(-1000, 1000, "e_" + std::to_string(i))) ;
        error_variables.push_back(solver.MakeNumVar(0, 1000, "l_" + std::to_string(i))) ;
        int n = error_variables.size() ;
        // error absolute value
        lp::MPConstraint* const abs_1 = solver.MakeRowConstraint(0.0, lp::MPSolver::infinity(), "abs_1_" + std::to_string(i));
        abs_1->SetCoefficient(error_variables[n - 1], 1) ;
        abs_1->SetCoefficient(error_variables[n - 2], +1.0) ;
        lp::MPConstraint* const abs_2 = solver.MakeRowConstraint(0.0, lp::MPSolver::infinity(), "abs_2_" + std::to_string(i));
        abs_2->SetCoefficient(error_variables[n - 1], 1) ;
        abs_2->SetCoefficient(error_variables[n - 2], -1.0) ;
        objective->SetCoefficient(error_variables[n - 1], 1) ;
        // kmer constraint
        int rhs = 50 + (i - 500 % 10) ;
        lp::MPConstraint* const ct = solver.MakeRowConstraint(rhs, rhs, "count_" + std::to_string(i));
        ct->SetCoefficient(track_variables[i % 10], 50) ;
        ct->SetCoefficient(error_variables[n - 2], 1.0) ;
    }
    objective->SetMinimization() ;
    solver.Solve() ;
}

void Genotyper::solve_lp() {
    cout << "Solving LP.." << endl ;
    auto c = Configuration::getInstance() ;
    int n = 0 ;
    int batch_size = 1 ;
    #pragma omp parallel for num_threads(48)
    for (int i = 0; i < clusters.size() / batch_size; i++) {
        genotype_clusters(i * batch_size, (i + 1) * batch_size) ;
        n += 1 ;
        if (n % 10 == 0) {
            cout << "Solved " << std::left << std::setw(6) << n << " out of " << std::setw(6) << clusters.size() / batch_size << " batches.\r" ;
        }
    }
    cout << "\n" ;
}

void Genotyper::genotype_clusters(int batch_b, int batch_e) {
    auto c = Configuration::getInstance() ;
    lp::MPSolver solver("cluster_" + std::to_string(batch_b), lp::MPSolver::CLP_LINEAR_PROGRAMMING) ;
    lp::MPObjective* const objective = solver.MutableObjective() ;
    vector<lp::MPVariable*> track_variables ;
    unordered_map<SimpleTrack, int> track_index ;
    int m = 0 ;
    for (int i = batch_b; i < batch_e; i++) {
        for (auto track = clusters[i].begin(); track != clusters[i].end(); track++) {
            if (tracks.find(*track) != tracks.end()) {
                track_variables.push_back(solver.MakeNumVar(0.0, 1.0, track->get_name())) ;
                track_index[*track] = track_variables.size() - 1 ;
                m += tracks[*track].size() ;
            }
        }
    }
    //cout << m << " kmers." << endl ;
    //cout << track_variables.size() << " track variables." << endl ;
    vector<lp::MPVariable*> error_variables ;
    unordered_map<uint64_t, bool> lp_kmers ;
    std::ofstream j ;
    //j.open(c->workdir + "/batch_" + std::to_string(batch_b) + ".json") ;
    string name = clusters[batch_b][0].get_name() ;
    j.open(c->workdir + "/batch_" + name + ".json") ;
    j << "{\n" ;
    int u = 0 ;
    //cout << "building lp.." << endl ;
    for (int i = batch_b; i < batch_e; i++) {
        for (auto track = clusters[i].begin(); track != clusters[i].end(); track++) {
            for (auto k = tracks[*track].begin(); k != tracks[*track].end(); k++) {
                if (lp_kmers.find(*k) == lp_kmers.end()) {
                    auto kmer = kmers[*k] ;
                    //print_kmer(kmer) ;
                    // JSON output
                    if (u != 0) {
                        j << ",\n" ;
                    }
                    j << "\"" << decode_kmer(*k) << "\":" ;
                    j << nlohmann::json(kmer).dump(4) ;
                    u += 1 ;
                    // JSON output
                    lp_kmers[*k] = true ;
                    int count = kmer.inverse ? kmer.total - kmer.count : kmer.count ;
                    int residue = kmer.inverse ? 0 : kmer.loci.size() - kmer.tracks.size() ;
                    int coverage = kmer.type == KMER_TYPE_JUNCTION ? this->coverage : gc_coverage[kmer.gc] ;
                    int lp_count = min(count, coverage * int(kmer.loci.size())) ;
                    //cout << kmer.count << " " << coverage << " " << lp_count << endl ;
                    error_variables.push_back(solver.MakeNumVar(-1000, 1000, "e_" + decode_kmer(*k))) ;
                    error_variables.push_back(solver.MakeNumVar(0, 1000, "l_" + decode_kmer(*k))) ;
                    int n = error_variables.size() ;
                    // error absolute value
                    lp::MPConstraint* const abs_1 = solver.MakeRowConstraint(0.0, lp::MPSolver::infinity(), "abs_1_" + decode_kmer(*k));
                    abs_1->SetCoefficient(error_variables[n - 1], 1) ;
                    abs_1->SetCoefficient(error_variables[n - 2], +1.0) ;
                    lp::MPConstraint* const abs_2 = solver.MakeRowConstraint(0.0, lp::MPSolver::infinity(), "abs_2_" + decode_kmer(*k));
                    abs_2->SetCoefficient(error_variables[n - 1], 1) ;
                    abs_2->SetCoefficient(error_variables[n - 2], -1.0) ;
                    objective->SetCoefficient(error_variables[n - 1], 1) ;
                    // kmer constraint
                    if (kmer.type == KMER_TYPE_JUNCTION || kmer.trend) {
                        int rhs = lp_count - residue * coverage ;
                        lp::MPConstraint* const ct = solver.MakeRowConstraint(rhs, rhs, "count_" + decode_kmer(*k));
                        for (auto _track = kmer.tracks.begin(); _track != kmer.tracks.end(); _track++) {
                            ct->SetCoefficient(track_variables[track_index[_track->first]], coverage * kmer.tracks[_track->first]) ;
                            ct->SetCoefficient(error_variables[n - 2], 1.0) ;
                        }
                    } else {
                        int rhs = lp_count - residue * coverage ;
                        for (auto _track = kmer.tracks.begin(); _track != kmer.tracks.end(); _track++) {
                            rhs -= coverage * kmer.tracks[_track->first] ;
                        }
                        lp::MPConstraint* const ct = solver.MakeRowConstraint(rhs, rhs, "count_" + decode_kmer(*k));
                        for (auto _track = kmer.tracks.begin(); _track != kmer.tracks.end(); _track++) {
                            ct->SetCoefficient(track_variables[track_index[_track->first]], -1 * coverage * kmer.tracks[_track->first]) ;
                            ct->SetCoefficient(error_variables[n - 2], 1.0) ;
                        }
                    }
                }
            }
        }
    }
    j << "\n}\n" ;
    //cout << lp_kmers.size() << " lp kmers." << endl ;
    //cout << error_variables.size() << " error variables." << endl ;
    objective->SetMinimization() ;
    // Export program
    std::ofstream o ;
    o.open(c->workdir + "/batch_" + name + ".lp") ;
    string program ;
    solver.ExportModelAsLpFormat(false, &program) ;
    o << program ; 
    solver.Solve() ;
    // load solution
    for (auto track = track_index.begin(); track != track_index.end(); track++) {
        genotyping_tracks[track->first].lp_value = track_variables[track->second]->solution_value() ;
    }
}

void Genotyper::export_genotypes() {
    cout << "Exporting genotypes.." << endl ;
    auto c = Configuration::getInstance() ;
    std::ofstream o ;
    o.open(c->workdir + "/genotypes.bed") ;
    o << "#CHROM\tBEGIN\tEND\tSVTPYE\tLP\tGENOTYPE\n" ;
    for (auto it = genotyping_tracks.begin(); it != genotyping_tracks.end(); it++) {
        auto t = it->first ;
        auto g = it->second ;
        string genotype ;
        if (g.lp_value > 0.75) {
            genotype = "1/1" ;
        } else if (g.lp_value > 0.25) {
            genotype = "1/0" ;
        } else {
            genotype = "0/0" ;
        }
        o << get_chromosome_name(t.chrom) << "\t" << t.begin << "\t" << t.end << "\t" << get_svtype(t.svtype) << "\t" << g.lp_value << "\t" << genotype << endl ;
    } 
}
