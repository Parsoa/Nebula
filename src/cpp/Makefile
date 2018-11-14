CC = g++
CFLAGS = 

encoder.o: encoder.cpp encoder.hpp
	$(CC) $(CFLAGS) -c encoder.cpp

main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp

main: main.o encoder.o
	$(CC) $(CFLAGS) -o a.out main.o encoder.o

counter: counter.cpp
	/software/gcc/8.2.0/lssc0-linux/bin/g++ -std=c++11 -O3 -o counter.out counter.cpp
string: counter_s.cpp
	/software/gcc/8.2.0/lssc0-linux/bin/g++ -std=c++11 -o counter_s.out counter_s.cpp
debug: counter.cpp
	/software/gcc/8.2.0/lssc0-linux/bin/g++ -std=c++11 -ggdb -o counter.out counter.cpp
