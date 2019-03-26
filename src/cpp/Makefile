CC = /software/gcc/8.2.0/lssc0-linux/bin/g++
CFLAGS = -O3 -funroll-loops -std=c++11 -I htslib/htslib -std=c++11
LDFLAGS = /software/htslib/1.8/lssc0-linux/lib/libhts.a -lz -lm -lpthread -llzma -lbz2 -lcurl -lcrypto
SOURCES = counter.cpp
OBJECTS = $(SOURCES:.cpp=.o)

counter.o: counter.cpp
	$(CC) -g $(CFLAGS) -c $^ -o $@

counter: counter.o
	$(CC) counter.o -o counter.out $(LDFLAGS)

debug: counter.cpp
	$(CC) -g $(CFLAGS) counter.cpp -o counter.out $(LDFLAGS)
