OR_TOOLS_DIR = or-tools
OR_TOOLS_LIB_DIR = $(OR_TOOLS_DIR)/lib
OR_TOOLS_INC_DIR = $(OR_TOOLS_DIR)/include
OR_TOOLS_FLAGS = -L$(shell pwd)/$(OR_TOOLS_LIB_DIR) -lglog -lprotobuf -lgflags -lCbc -lCbcSolver -lClp -lClpSolver -lCoinUtils -lOsiCbc -lOsiClp -lCgl -lortools

HTSLIB_PREFIEX = ./htslib
HTSLIBFLAGS = -L$(HTSLIB_PREFIEX) -lz -llzma -lbz2 -lcurl -lcrypto -lhts
HTSLIBSOURCES = htslib

CXX = g++
LDFLAGS = -lm -lpthread $(HTSLIBFLAGS)
CXXFLAGS = -std=c++11 -O3 -funroll-loops -fopenmp -I$(HTSLIB_PREFIEX)

HPP = $(wildcard *.hpp)
SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)

LINEARSOL ?= 0
ifneq (${LINEARSOL}, 0)
	CXXFLAGS += -I $(OR_TOOLS_INC_DIR)
	LDFLAGS += $(OR_TOOLS_FLAGS)
	OR_TOOLS_TARGET ?= .ortools
endif

.ortools:
	wget https://github.com/google/or-tools/releases/download/v7.8/or-tools_ubuntu-16.04_v7.8.7959.tar.gz; mkdir or-tools; tar -xzvf or-tools_ubuntu-16.04_v7.8.7959.tar.gz -C or-tools --strip-components=1

debug: CXXFLAGS += -DDEBUG -g

all: nebula
debug: nebula
docker: nebula

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

nebula: $(OBJ) $(HPP) $(OR_TOOLS_TARGET)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(LDFLAGS)

clean:
	rm *.o
	rm -rf or-tools*
.PHONY: clean

