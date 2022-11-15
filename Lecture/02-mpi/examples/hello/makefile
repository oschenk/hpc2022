all: hello hello_advanced hello_critical

hello: hello.c
	mpicc -O3 $< -o $@

hello_advanced: hello_advanced.c
	mpicc -O3 $< -o $@

hello_critical: hello_critical.c
	mpicc -O3 $< -o $@

clean:
	rm -rf hello hello_advanced hello_critical

