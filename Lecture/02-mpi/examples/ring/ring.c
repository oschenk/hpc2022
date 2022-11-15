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

	int right = (my_rank + 1)        % size;
	int left  = (my_rank - 1 + size) % size;

	/* ... this SPMD-style neighbor computation with modulo has the same meaning as: */
	/* right = my_rank + 1;          */
	/* if (right == size) right = 0; */
	/* left = my_rank - 1;           */
	/* if (left == -1) left = size-1;*/

	int sum = 0;
	for (int i = 0, sent_buf = my_rank, recv_buf; i < size; i++, sent_buf = recv_buf)
	{
		MPI_Request request;
		MPI_ERR_CHECK(MPI_Issend(&sent_buf, 1, MPI_INT, right, to_right,
			MPI_COMM_WORLD, &request));

		MPI_Status status;
		MPI_ERR_CHECK(MPI_Recv(&recv_buf, 1, MPI_INT, left, to_right,
			MPI_COMM_WORLD, &status));

		MPI_ERR_CHECK(MPI_Wait(&request, &status));

		sum += recv_buf;
	}

	printf ("PE%i:\tSum = %i\n", my_rank, sum);

	MPI_ERR_CHECK(MPI_Finalize());

	return 0;
}

