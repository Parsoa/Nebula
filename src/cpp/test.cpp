#include <locale>
#include <ctime>
#include <string>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <cstdint>
#include <string>
#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>

void handler(int sig) {
    void *array[10];
    size_t size;
    // get void*'s for all entries on the stack
    size = backtrace(array, 10);
    // print out all the frames to stderr
    fprintf(stderr, "Error: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    exit(1);
}
int main(int argc, char** argv) {
    //signal(SIGSEGV, handler);
    signal(SIGSEGV, handler);
    std::unordered_map<std::string, std::vector<int>> index ;
    index.insert(std::make_pair("AAAAAAAAAAAAAAC", std::vector<int>(1, 2))) ;
    index.insert(std::make_pair("CAAAAAAAAAAAAAA", std::vector<int>(1, 2))) ;
}
