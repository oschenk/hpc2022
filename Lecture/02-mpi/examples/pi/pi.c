#include <limits.h>
#include <mpi.h>
#include <stdio.h>
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
    // Start MPI application.
    double pi;
    MPI_ERR_CHECK(MPI_Init(&argc, &argv));

    if (argc != 2)
    {
        printf("Usage: %s <n>\n", argv[0]);
        exit(1);
    }

    // Get the N argument value.	
    size_t n = (size_t)strtoul(argv[1], NULL, 0);

    // Start timer.
    double start = MPI_Wtime();

    // Get current process rank and number of processes in
    // MPI application.
    int rank, size;
    MPI_ERR_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
    MPI_ERR_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &size));

    // Compute the number of points per process.
    size_t npoints = n / size;
    size_t npointsr = n / size;
    if (rank == size - 1)
        npointsr += n % size;

    // Calculate process' part of intergral.
    double w = 1.0 / n;
    double psum = 0.0;	
    for (size_t i = rank * npoints, j = 0; j < npointsr; i++, j++)
    {
        double x = w * ((double)i + 0.5);
        psum += 4.0 / (1.0 + x * x);
    }

    // Sum partial integrals from MPI processes
    if (rank == 0) 
    {
        double sum = 0.0;
        MPI_ERR_CHECK(MPI_Reduce(&psum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD));
        pi = w * sum;
    }
    else
    {
        MPI_ERR_CHECK(MPI_Reduce(&psum, NULL, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD));
    }


    // Stop timer.
    double finish = MPI_Wtime();

    // Print results on rank 0.
    if (rank == 0)
        printf("Calculated pi = %24.16g on %d MPI processors in %f secs\n", pi, size, finish - start);

    // Finish MPI application.
    MPI_ERR_CHECK(MPI_Finalize());

    return 0;
}

