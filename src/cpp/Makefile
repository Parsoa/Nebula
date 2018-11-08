CC = g++
CFLAGS = 

encoder.o: encoder.cpp encoder.hpp
	$(CC) $(CFLAGS) -c encoder.cpp

main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp

main: main.o encoder.o
	$(CC) $(CFLAGS) -o a.out main.o encoder.o

# all: c_counter.cpp
# 	/software/gcc/8.2.0/lssc0-linux/bin/g++ -std=c++11 -o c_counter.out c_counter.cpp
# debug: c_counter.cpp
# 	/software/gcc/8.2.0/lssc0-linux/bin/g++ -std=c++11 -g -o c_counter.out c_counter.cpp
