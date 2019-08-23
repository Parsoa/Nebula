CC = g++
CFLAGS = -O3 -funroll-loops -std=c++11 -I htslib
LDFLAGS = /software/htslib/1.8/lssc0-linux/lib/libhts.a -lz -lm -lpthread -llzma -lbz2 -lcurl -lcrypto
# Docker config
DCFLAGS = -O3 -funroll-loops -std=c++11
DLDFLAGS = -lz -lm -lpthread -lhts
SOURCES = counter.cpp
OBJECTS = $(SOURCES:.cpp=.o)

counter.o: counter.cpp
	$(CC) -g $(CFLAGS) -c $^ -o $@

counter: counter.o
	$(CC) counter.o -o counter.out $(LDFLAGS)

debug: counter.cpp
	$(CC) -g $(CFLAGS) -DDEBUG counter.cpp -o counter.out $(LDFLAGS)

docker: counter.cpp
	$(CC) -g $(DCFLAGS) counter.cpp -o counter.out $(DLDFLAGS)

