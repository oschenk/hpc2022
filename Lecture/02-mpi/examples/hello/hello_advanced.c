#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <unistd.h>


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

    int my_rank, size;
    char hostname[256];
    MPI_ERR_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &my_rank));
    MPI_ERR_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &size));
    gethostname(hostname,255);

    if (my_rank == 0)   
    {
        printf ("I am process %i out of %i on host %s: Hello world!\n",
                my_rank, size, hostname);
    }
    else
    {
        printf("I am process %i out of %i on host %s\n", my_rank, size,
                hostname);
    }

    MPI_ERR_CHECK(MPI_Finalize());

    return 0;
}

