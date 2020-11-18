OR_TOOLS_DIR = /share/hormozdiarilab/Codes/NebulousSerendipity/binaries/or-tools_Ubuntu-16.04-64bit_v7.8.7959
docker: OR_TOOLS_DIR = /lib/or-tools_Ubuntu-16.04-64bit_v7.8.7959
OR_TOOLS_LIB_DIR = $(OR_TOOLS_DIR)/lib
OR_TOOLS_INC_DIR = $(OR_TOOLS_DIR)/include
OR_TOOLS_FLAGS = -L$(OR_TOOLS_LIB_DIR) -lglog -lprotobuf -lgflags -lCbc -lCbcSolver -lClp -lClpSolver -lCoinUtils -lOsiCbc -lOsiClp -lCgl -lortools

HTSLIBFLAGS = /software/htslib/1.8/lssc0-linux/lib/libhts.a -lz -llzma -lbz2 -lcurl -lcrypto
docker: HTSLIBFLAGS = -lz -llzma -lbz2 -lhts

CXX = g++
LDFLAGS = -lm -lpthread $(HTSLIBFLAGS) $(OR_TOOLS_FLAGS)
CXXFLAGS = -std=c++11 -O3 -funroll-loops -fopenmp -I htslib -I $(OR_TOOLS_INC_DIR)

HPP = $(wildcard *.hpp)
SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)

all: nebula

debug: CXXFLAGS += -DDEBUG -g

docker: nebula

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

nebula: $(OBJ) $(HPP) 
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(LDFLAGS)

debug: $(OBJ) $(HPP) 
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(LDFLAGS)

clean:
	rm *.o
.PHONY: clean

