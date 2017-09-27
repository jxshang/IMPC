all: IMPC

IMPC : IMPC.o Graph.o Node.o
	g++ -O3 -o IMPC IMPC.o Graph.o Node.o

IMPC.o : IMPC.cpp
	g++ -O3 -c IMPC.cpp

Graph.o : Graph.cpp
	g++ -O3 -c Graph.cpp

Node.o : Node.cpp
	g++ -O3 -c Node.cpp

clean :
	rm Node.o Graph.o IMPC.o IMPC
