CC = g++
CFLAGS = -std=c++11 -O3 -funroll-loops -I htslib -fopenmp
LDFLAGS = /software/htslib/1.8/lssc0-linux/lib/libhts.a -lz -lm -lpthread -llzma -lbz2 -lcurl -lcrypto
# Docker config
DCFLAGS = -O3 -funroll-loops -std=c++11
DLDFLAGS = -lz -lm -lpthread -lhts
SOURCES = counter.cpp
OBJECTS = $(SOURCES:.cpp=.o)

debug: counter.cpp
	$(CC) -g $(CFLAGS) -DDEBUG_MODE $< -o counter.out $(LDFLAGS)

docker: counter.cpp
	$(CC) $(DCFLAGS) $< -o counter.out $(DLDFLAGS)

counter: counter.cpp
	$(CC) -g $(CFLAGS) $< -o counter.out $(LDFLAGS)
