kmc.so: kmc.o
	g++ -shared -o kmc.so kmc.o

kmc.out: kmc.o
	g++ *.o -o kmc.out

kmc.o: kmc.cpp
	g++ -c -std=c++14 -Wall -fpic *.cpp

clean:
	rm *.so *.o *.out
