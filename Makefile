all: IMPC

IMPC : IMPC.o
	g++ -O3 -o IMPC IMPC.o

IMPC.o : IMPC.cpp
	g++ -O3 -c IMPC.cpp

clean :
	rm IMPC.o IMPC
