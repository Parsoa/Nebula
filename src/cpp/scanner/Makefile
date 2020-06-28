CC = g++
CFLAGS = -std=c++11 -O3 -funroll-loops -fopenmp
LDFLAGS = -lm -lpthread -lcurl
SOURCES = scanner.cpp
OBJECTS = $(SOURCES:.cpp=.o)
# Docker config
DCFLAGS = -O3 -funroll-loops -std=c++11
DLDFLAGS = -lm -lpthread

debug: scanner.cpp
	$(CC) -g $(CFLAGS) -DDEBUG_MODE $< -o scanner.out $(LDFLAGS)

shared: scanner.cpp
	$(CC) $(CFLAGS) -fPIC -shared $< -o libnebula.so $(LDFLAGS)

docker: scanner.cpp
	$(CC) $(DCFLAGS) $< -o scanner.out $(DLDFLAGS)

scanner: scanner.cpp
	$(CC) -g $(CFLAGS) $< -o scanner.out $(LDFLAGS)
