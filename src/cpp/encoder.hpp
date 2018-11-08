#include <string>
#include <cstdint>

class ReadEncoder {

public:

    uint8_t* encode(std::string read) ;

    void decode() ;

private:

    uint8_t* bead ;

};