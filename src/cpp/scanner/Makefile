CXX = g++
LDFLAGS = -lm -lpthread -lcurl
CXXFLAGS = -std=c++11 -O3 -funroll-loops -fopenmp

HPP = $(wildcard *.hpp)
SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)

all: scanner 

debug: CXXFLAGS += -DDEBUG -g
debug: preprocessor

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

scanner: $(OBJ) $(HPP) 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJ)

clean:
	rm *.o
.PHONY: clean

