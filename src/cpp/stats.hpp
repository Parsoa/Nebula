#ifndef STT_HPP
#define STT_HPP

class NormalDistribution {

public:

    NormalDistribution(double mean, double std): mean(mean), std(std), var(std * std) {}
    double prob(double x) ;
    double log_prob(double x) ;

private:

    double var ;
    double std ;
    double mean ;

} ;

#endif
