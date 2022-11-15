all: ring 

ring: ring.c
	mpicc -std=c99 -O3 $< -o $@

clean:
	rm -rf ring 

