all: ring ring_advanced_irecv_issend ring_advanced_irecv_ssend

ring: ring.c
	mpicc -std=c99 -O3 $< -o $@

ring_advanced_irecv_issend: ring_advanced_irecv_issend.c
	mpicc -std=c99 -O3 $< -o $@

ring_advanced_irecv_ssend: ring_advanced_irecv_ssend.c
	mpicc -std=c99 -O3 $< -o $@

clean:
	rm -rf ring ring_advanced_irecv_issend ring_advanced_irecv_ssend

