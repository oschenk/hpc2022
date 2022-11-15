/****************************************************************
 *                                                              *
 * This file has been written as a sample solution to an        *
 * exercise in a course given at the High Performance           *
 * Computing Centre Stuttgart (HLRS).                           *
 * The examples are based on the examples in the MPI course of  *
 * the Edinburgh Parallel Computing Centre (EPCC).              *
 * It is made freely available with the understanding that      *
 * every copy of this file must include this header and that    *
 * HLRS and EPCC take no responsibility for the use of the      *
 * enclosed teaching material.                                  *
 *                                                              *
 * Authors: Joel Malard, Alan Simpson,            (EPCC)        *
 *          Rolf Rabenseifner, Traugott Streicher (HLRS)        *
 *                                                              *
 * Contact: rabenseifner@hlrs.de                                * 
 *                                                              *  
 * Purpose: A program to try MPI_Issend and MPI_Recv.           *
 *                                                              *
 * Contents: C-Source                                           *
 *                                                              *
 ****************************************************************/

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

#define MPI_ERR_CHECK(call) {                        \
    do { int err = call; if (err != MPI_SUCCESS) {   \
        char errstr[MPI_MAX_ERROR_STRING];           \
        int szerrstr;                                \
        MPI_Error_string(err, errstr, &szerrstr);    \
        fprintf(stderr, "MPI error at %s:%i : %s\n", \
            __FILE__, __LINE__, errstr);             \
        MPI_Abort(MPI_COMM_WORLD, 1);                \
    }} while (0);                                    \
}

#define to_right 201

int main(int argc, char* argv[])
{
	MPI_ERR_CHECK(MPI_Init(&argc, &argv));

	int my_rank, size;
	MPI_ERR_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &my_rank));
	MPI_ERR_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &size));

        // TODO Compute rank of neighboring processes
	int right = 0;
	int left  = 0;

	int sum = 0;
	for (int i = 0, sent_buf = my_rank, recv_buf; i < size; i++, sent_buf = recv_buf)
	{
		MPI_Request request;
		MPI_Status status;

                // TODO Send sent_buf to your right neighbor
                // TODO Get buffer from your left neighbor into recv_buf

		sum += recv_buf;
	}

	printf ("PE%i:\tSum = %i\n", my_rank, sum);

	MPI_ERR_CHECK(MPI_Finalize());

	return 0;
}

