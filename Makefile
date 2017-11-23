coloring: ./src/coloring.cpp
	gcc ./src/graphio.c -c -O3
	gcc ./src/mmio.c -c -O3
#	g++ ./src/GraphColoring.cpp -c -O3 -std=c++14
	g++ ./src/coloring.cpp -fopenmp -c -O3 -std=c++14
#	g++ -o coloring coloring.o mmio.o graphio.o GraphColoring.o -O3 -fopenmp -std=c++14
	g++ -o coloring coloring.o mmio.o graphio.o -O3 -fopenmp -std=c++14
clean:
	rm coloring *.o
