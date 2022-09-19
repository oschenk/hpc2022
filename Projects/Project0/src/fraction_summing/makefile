# by default select gcc
CXX=g++
CXXFLAGS=-O3


SOURCES = fraction_toolbox.cpp
HEADERS = fraction_toolbox.hpp
OBJ     = fraction_toolbox.o

.SUFFIXES: .cpp

all: main

fraction_toolbox.o: fraction_toolbox.cpp fraction_toolbox.hpp
	$(CXX) $(CXXFLAGS) -c fraction_toolbox.cpp -o fraction_toolbox.o

main: $(OBJ) main.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) *.o main.cpp -o main

clean:
	rm -f main
	rm -f *.o
	rm -f *.i
	rm -f *.lst
	rm -f output.*