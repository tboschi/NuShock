d1:	mcflux.cxx ../eventrate/flux.h ../eventrate/flux.c
	g++ -g -std=c++11 -c mctools.c -o mctools.o -I. -I../eventrate
	g++ -g -std=c++11 -c mcflux.cxx -o mcflux.o -I. -I../eventrate
	g++ -g -std=c++11 -c ../eventrate/flux.c -o flux.o -I. -I../eventrate
	g++ -g -std=c++11 -o mcflux mctools.o mcflux.o flux.o -lgomp -lnlopt -lgsl -lgslcblas
	rm *.o

