all: pingpong

pingpong: pingpong.c
	mpicc -std=c99 -O3 $< -o $@

clean:
	rm -rf pingpong

