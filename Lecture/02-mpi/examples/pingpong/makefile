all: pingpong pingpong_advanced_send pingpong_advanced_ssend

pingpong: pingpong.c
	mpicc -std=c99 -O3 $< -o $@

pingpong_advanced_send: pingpong_advanced_send.cpp
	mpiCC -O3 $< -o $@

pingpong_advanced_ssend: pingpong_advanced_ssend.cpp
	mpiCC -O3 $< -o $@

clean:
	rm -rf pingpong pingpong_advanced_send pingpong_advanced_ssend

