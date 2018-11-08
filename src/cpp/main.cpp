#include "encoder.hpp"

int main() {
    ReadEncoder* r = new ReadEncoder() ;
    r->encode("GGATAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACAGAGATAACTATTGATACAACACCTTCATGACCC") ;
    r->decode() ;
}