#include <omp.h>
#include <errno.h>
#include <iomanip>

#include "preprocessor.hpp"

#include "inner.hpp"
#include "logger.hpp"
#include "junction.hpp"
#include "bed_utils.hpp"

using namespace std ;

void Preprocessor::run() {
    auto c = Configuration::getInstance() ;
    cout << "Preprocessing.." << endl ;
    load_chromosomes(c->reference) ;
    auto inner_kmer_extractor = InnerKmerExtractor(c->threads) ;
    inner_kmer_extractor.run() ;
    auto junction_kmer_extractor = JunctionKmerExtractor(c->threads) ;
    junction_kmer_extractor.run() ;
    merge_genotyping_kmers(inner_kmer_extractor.inner_kmers, junction_kmer_extractor.junction_kmers) ;
    scan_reference(c->threads) ;
    filter_kmers() ;
    dump_kmers(c->workdir) ;
}

// ============================================================================= \\
// ============================================================================= \\
// ================================ FASTA Files ================================ \\
// ============================================================================= \\
// ============================================================================= \\

void Preprocessor::scan_reference(int threads) {
    cout << "--------------------------------------------------------- " << endl ;
    cout << "Scanning reference genome.." << endl ;
    threads = min(threads, int(chromosome_seqs.size())) ;
    for (int t = 0; t < threads - 1; t++) {
        cout << endl ;
    }
    int m = 0 ;
    int n = chromosomes.size() ;
    vector<unordered_map<uint64_t, Kmer>> batches(chromosomes.size()) ; 
    while (m < n) {
        int p = m ;
        #pragma omp parallel for num_threads(threads)
        for (int i = 0; i < threads; i++) {
            if (m + i < n) {
                batches[i + m] = scan_chromosome(chromosomes[m + i], threads) ;
            }
        }
        m += threads ;
        if (m >= chromosomes.size()) {
            m = n ;
        }
    }
    cout << "--------------------------------------------------------- " << endl ;
    for (int i = 0; i < chromosomes.size(); i++) {
        cout << "Merging " << std::left << std::setw(6) << chromosomes[i] << ", " << std::setw(9) << batches[i].size() << " matches.." << endl ;
        for(auto it = batches[i].begin(); it != batches[i].end(); it++) {
            kmers[it->first].count += it->second.count ;
            kmers[it->first].loci.insert(kmers[it->first].loci.end(), it->second.loci.begin(), it->second.loci.end()) ;
        }
        batches[i].clear() ;
    }
    cout << "Reference scan completed." << endl ;
}

unordered_map<uint64_t, Kmer> Preprocessor::scan_chromosome(string chrom, int threads) {
    unordered_map<uint64_t, Kmer> _kmers ;
    uint64_t k = 0 ;
    uint64_t left = 0 ;
    uint64_t right = 0 ;
    uint64_t l = strlen(chromosome_seqs[chrom]) ;
    //cout << "Scanning " << chrom << " with " << l << " bases.." << endl ;
    KmerIterator it(chromosome_seqs[chrom], 234, l - 500, 234, ITER_MODE_REF) ;
    int n = 0 ;
    while (it) {
        auto k = it->kmer ;
        if (kmers.find(k) != kmers.end()) {
            if (_kmers.find(k) == _kmers.end()) {
                _kmers[k] = Kmer() ;
                _kmers[k].count = 0 ;
            }
            if (_kmers[k].count >= 4) {
                it++ ;
                continue ;
            }
            Locus locus({get_chromosome_index(chrom), uint32_t(it.position), LOCUS_TYPE_REF, it->left, it->right, it->gc / 5, false}) ;
            _kmers[k].loci.push_back(locus) ;
            _kmers[k].count ++ ;
            n++ ;
        } else { 
            uint64_t rc = encoded_reverse_complement(k) ;
            if (kmers.find(rc) != kmers.end()) {
                if (_kmers.find(rc) == _kmers.end()) {
                    _kmers[rc] = Kmer() ;
                    _kmers[rc].count = 0 ;
                }
                if (_kmers[rc].count >= 4) {
                    it++ ;
                    continue ;
                }
                Locus locus({get_chromosome_index(chrom), uint32_t(it.position), LOCUS_TYPE_REF, it->left, it->right, it->gc / 5, false}) ;
                _kmers[rc].loci.push_back(locus) ;
                _kmers[rc].count ++ ;
                n++ ;
            }
        }
        it++ ;
        if (it.position % 10000000 == 0) {
            cout_mutex.lock() ;
            int index = get_chromosome_index(chrom) ;
            for (int j = 0; j < (threads - 1) - index; j++) {
                cout << "\x1b[A" ;
            }
            cout << "\r"<< std::left << std::setw(6) << chrom << " progress " << std::fixed << std::setprecision(3) << float(it.position) / float(l) << "%, loci: " << std::left << std::setw(9) << n ;
            for (int j = 0; j < (threads - 1) - index; j++) {
                cout << endl ;
            }
            cout_mutex.unlock() ;
            //cout << chrom << " progress " << float(it.position) / float(l) << "%, loci: " << n << endl ;
        }
    }
    return _kmers ;
}

// ============================================================================= \\
// ================================= Filtering ================================= \\
// ============================================================================= \\

void Preprocessor::merge_genotyping_kmers(unordered_map<uint64_t, Kmer> inner_kmers, unordered_map<uint64_t, Kmer> junction_kmers) {
    auto c = Configuration::getInstance() ;
    cout << "--------------------------------------------------------- " << endl ;
    cout << "Merging inner and junction kmers.." << endl ;
    std::vector<uint64_t> remove ;
    for (auto junction_kmer = junction_kmers.begin(); junction_kmer != junction_kmers.end(); junction_kmer++) {
        if (inner_kmers.find(junction_kmer->first) != inner_kmers.end()) {
            remove.push_back(junction_kmer->first) ;
        }
    }
    int m = 0 ;
    std::ofstream j ;
    j.open(c->workdir + "/common.json") ;
    j << "{\n" ;
    for (auto kmer = remove.begin(); kmer != remove.end(); kmer++) {
        if (m != 0) {
            j << ",\n" ;
        }
        m++ ;
        j << "\"" << decode_kmer(*kmer) << "\":" ;
        j << nlohmann::json(inner_kmers.find(*kmer)->second).dump(4) ;
        inner_kmers.erase(*kmer) ;
        junction_kmers.erase(*kmer) ;
    }
    j << "\n}\n" ;
    j.close() ;
    cout << "Filtered " << m << " common kmers." << endl ;
    kmers.insert(inner_kmers.begin(), inner_kmers.end()) ;
    kmers.insert(junction_kmers.begin(), junction_kmers.end()) ;
}

// ============================================================================= \\
// ================================= Filtering ================================= \\
// ============================================================================= \\

vector<uint64_t> find_interest_masks(Kmer& k) {
    std::vector<uint64_t> interset_masks ;
    for (auto locus = k.loci.begin(); locus != k.loci.end(); locus++) {
        if (locus->type != LOCUS_TYPE_REF) { // these will be junction
            if (locus->left != 0) {
                interset_masks.push_back(locus->left) ;
            }
            if (locus->right != 0) {
                interset_masks.push_back(locus->right) ;
            }
        }
        // these will be inner and breakpoint loci from deletions + inner breakpoint from insertions
        else {
            for (auto track = k.tracks.begin(); track != k.tracks.end(); track++) {
                if (locus->chrom == track->first.chrom) {
                    if (locus->position >= track->first.begin - 32 && locus->position <= track->first.end) {
                        if (locus->left != 0) {
                            interset_masks.push_back(locus->left) ;
                        }
                        if (locus->right != 0) {
                            interset_masks.push_back(locus->right) ;
                        }
                        locus->type = LOCUS_TYPE_INNER ;
                    }
                }
            }
        }
    }
    bool found = false ;
    for (auto locus = k.loci.begin(); locus != k.loci.end(); locus++) {
        if (locus->type != LOCUS_TYPE_REF) {
            found = true ;
        }
    }
    if (not found) {
        print_kmer(k) ;
    }
    assert(found) ;
    return interset_masks ;
}

void Preprocessor::filter_kmers() {
    auto c = Configuration::getInstance() ;
    int l = 0 ;
    int n = 0 ;
    int r = 0 ;
    int p = 0 ;
    int threads = min(c->threads, int(kmers.bucket_count())) ;
    vector<vector<uint64_t>> filters(threads) ;
    cout << "Filtering " << kmers.size() << " kmers using " << threads << " threads on " << kmers.bucket_count() << " buckets.." << endl ;
    #pragma omp parallel for num_threads(threads)
    for (size_t bucket = 0; bucket < kmers.bucket_count(); bucket++) {
        //cout << bucket << " " << n << endl ;
        for (auto kmer = kmers.begin(bucket); kmer != kmers.end(bucket); kmer++) {
            n++ ;
            Kmer& k = kmer->second ;
            // filter based on ref count
            if (k.count > 3) {
                int t = omp_get_thread_num() ;
                filters[t].push_back(kmer->first) ;
                r += 1 ;
                //cout << "filtered based on ref count" << endl ;
                continue ;
            }
            std::vector<Locus> loci(k.loci) ;
            std::vector<uint64_t> interset_masks = find_interest_masks(k) ;
            bool found = false ;
            for (auto locus = k.loci.begin(); locus != k.loci.end(); locus++) {
                if (locus->type != LOCUS_TYPE_REF) {
                    found = true ;
                }
            }
            if (not found) {
                print_kmer(k) ;
            }
            assert(found) ;
            // number of loci before filtering
            int l_1 = k.loci.size() ;
            auto locus = k.loci.begin() ;
            while (locus != k.loci.end()) {
                if (locus->type == LOCUS_TYPE_REF) {
                    bool found = false ;
                    for (auto m = interset_masks.begin(); m != interset_masks.end(); m++) {
                        if (locus->left != 0) {
                            if (is_canonical_subsequence(locus->left, *m)) {
                                found = true ;
                                break ;
                            }
                        }
                        if (locus->right != 0) {
                            if (is_canonical_subsequence(locus->right, *m)) {
                                found = true ;
                                break ;
                            }
                        }
                    }
                    // doesn't have any shared masks, filter
                    if (not found) {
                        k.filtered_loci.push_back(*locus) ;
                        locus = k.loci.erase(locus) ;
                        continue ;
                    }
                }
                locus++ ;
            }
            // number of loci after filtering
            int l_2 = k.loci.size() ;
            //cout << l_2 << " remaining" << endl ;
            // loci with less than two masks exist
            if (k.type == KMER_TYPE_JUNCTION) {
                auto a = find_if(k.loci.begin(), k.loci.end(), [](Locus l) {
                    return l.left == 0 || l.right == 0 ;
                }) ;
                //non-junction loci exists
                auto b = find_if(k.loci.begin(), k.loci.end(), [](Locus l) {
                    return l.type == LOCUS_TYPE_REF ;
                }) ;
                // won't happen for inner kmers
                if (a != k.loci.end()) { // junction loci exist wtih less than two masks
                    if (l_1 != l_2) { // some ref loci were filtered
                        if (b == k.loci.end()) { // all ref loci were filtered
                            // because ref loci will have masks, count them instead and subtract from total. May overcount.
                            // Counting with one mask will undercount
                            k.loci = loci ;
                            auto locus = k.loci.begin() ; // count non-junction loci only
                            while (locus != k.loci.end()) {
                                if (locus->type == LOCUS_TYPE_JUNCTION) {
                                    k.junction_loci.push_back(*locus) ;
                                    locus = k.loci.erase(locus) ;
                                } else {
                                    locus++ ;
                                }
                            }
                            k.inverse = true ;
                        } else { // some ref loci remain
                            int t = omp_get_thread_num() ;
                            filters[t].push_back(kmer->first) ;
                            l += 1 ;
                            continue ;
                        }
                    } else { 
                        // no ref loci was filtered, so we need to count every loci
                    }
                }
            }
        }
        if (n - p > 10000) {
            p = n ;
            cout << "Progress " << std::left << std::setw(6) << std::fixed << std::setprecision(3) << float(n) / kmers.size() << "%..\r" << flush ;
        }
    }
    cout << endl ;
    cout << "Removing filtered kmers.." << endl ;
    int f = 0 ;
    std::ofstream j ;
    j.open(c->workdir + "/filtered.json") ;
    j << "{\n" ;
    //TODO: this is not necessary
    for (auto it = filters.begin(); it != filters.end(); it++) {
        for (auto k: *it) {
            if (f != 0) {
                j << ",\n" ;
            }
            f++ ;
            j << "\"" << decode_kmer(k) << "\":" ;
            j << nlohmann::json(kmers.find(k)->second).dump(4) ;
            kmers.erase(kmers.find(k)) ;
        }
    }
    cout << "Removed " << f << " kmers." << endl ;
    j << "\n}\n" ;
    j.close() ;
    cout << "Remaining " << kmers.size() << endl ;
}

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

void Preprocessor::dump_kmers(string path) {
    auto c = Configuration::getInstance() ;
    cout << "Verifying and dumping kmers.." << endl ;
    vector<int> counters ;
    int num_batches = 100 ;
    vector<mutex> locks(num_batches) ;
    vector<ofstream> output_files ;

    // this is for internal testing
    //ofstream test(path + "/kmers.txt") ;
    //for (auto kmer = kmers.begin(); kmer != kmers.end(); kmer++) {
    //    assert(kmer->second.weight == 1.0) ;
    //    test << decode_kmer(kmer->first) << endl ;
    //}
    //test.close() ;

    for (int i = 0; i < num_batches; i++) {
        string p = path + "/kmers_batch_" + std::to_string(i) + ".json" ;
        output_files.emplace_back(ofstream {p}) ;
        output_files[i] << "{\n" ;
        counters.push_back(0) ;
    }
    //TODO: check if correctly parallel
    //#pragma omp parallel for num_threads(c->threads)
    for (size_t bucket = 0; bucket < kmers.bucket_count(); bucket++) {
        int t = bucket % num_batches ;
        locks[t].lock() ;
        for (auto kmer = kmers.begin(bucket); kmer != kmers.end(bucket); kmer++) {
            bool found = false ;
            for (auto locus = kmer->second.loci.begin(); locus != kmer->second.loci.end(); locus++) {
                if (locus->type != LOCUS_TYPE_REF) {
                    found = true ;
                }
            }
            if (!found) {
                assert(kmer->second.inverse) ;
            }
            if (kmer->second.inverse) {
                assert(kmer->second.type == KMER_TYPE_JUNCTION) ;
            }
            if (counters[t] != 0) {
                output_files[t] << ",\n" ;
            }
            counters[t] += 1 ;
            output_files[t] << "\"" << decode_kmer(kmer->first) << "\":" ;
            output_files[t] << nlohmann::json(kmer->second).dump(4) ;
        }
        locks[t].unlock() ;
    }
    for (int i = 0; i < num_batches; i++) {
        output_files[i] << "\n}\n" ;
    }
}

