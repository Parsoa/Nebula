#include <omp.h>
#include <vector>
#include <mutex>
#include <iostream>
#include <chrono>
#include <pthread.h>
#include <thread>

using namespace std ;

void test_print() {
    int threads = 24 ;
    cout << "Testing with " << threads << " threads.." << endl ;
    std::mutex cout_mutex ;
    for (int t = 0; t < threads; t++) {
        cout << endl ;
    }
    cout << endl ;
    std::vector<uint64_t> counts(threads) ;
    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < threads; i++) {
        while (true) {
            int t = omp_get_thread_num() ;
            counts[t] += 1 ;
            if (counts[t] % (1000 + i) == 0) {
                cout_mutex.lock() ;
                for (int j = 0; j <= threads - t; j++) {
                    cout << "\x1b[A" ;
                }
                cout << "\r" ;
                cout << "Thread " << t << " printing.." ;
                for (int j = 0; j < threads - t; j++) {
                    cout << endl ;
                }
                cout_mutex.unlock() ;
            }
        }
    }
}

void test_print_2() {
    int threads = 4 ;
    cout << "Testing with " << threads << " threads.." << endl ;
    for (int t = 0; t < threads - 1; t++) {
        cout << t << endl ;
    }
    //cout << "\x1b[A" ;
    std::this_thread::sleep_for (std::chrono::seconds(2)) ;
    int t = 0 ;
    int n = 0 ;
    while (true) {
        for (int j = 0; j < (threads - 1) - t; j++) {
            cout << "\x1b[A" ;
        }
        //std::this_thread::sleep_for (std::chrono::seconds(2)) ;
        //cin.get() ;
        cout << "\rThread " << t << " printing " << n ;
        //cin.get() ;
        //std::this_thread::sleep_for (std::chrono::seconds(2)) ;
        for (int j = 0; j < (threads - 1) - t; j++) {
            cout << endl ;
        }
        t += 1 ;
        t %= threads ;
        n += 1 ;
    }
}
int main(int argc, char** argv) {
    for (int i = 100; i < 180; i++) {
        cout << i << endl ;
    }
    test_print_2() ;
}

