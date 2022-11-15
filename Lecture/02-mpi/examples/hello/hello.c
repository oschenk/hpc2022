#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

#define MPI_ERR_CHECK(call)                          \
    do { int err = call; if (err != MPI_SUCCESS) {   \
        char errstr[MPI_MAX_ERROR_STRING];           \
        int szerrstr;                                \
        MPI_Error_string(err, errstr, &szerrstr);    \
        fprintf(stderr, "MPI error at %s:%i : %s\n", \
            __FILE__, __LINE__, errstr);             \
        MPI_Abort(MPI_COMM_WORLD, 1);                \
    }} while (0);

int main(int argc, char* argv[])
{
	MPI_ERR_CHECK(MPI_Init(&argc, &argv));

	printf("Hello world!\n");

	MPI_ERR_CHECK(MPI_Finalize());
   
	return 0;
}

