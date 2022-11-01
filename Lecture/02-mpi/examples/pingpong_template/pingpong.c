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
 * Purpose: A program to try MPI_Ssend and MPI_Recv.            *
 *                                                              *
 * Contents: C-Source                                           *
 *                                                              *
 ****************************************************************/

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

#define proc_A 0
#define proc_B 1
#define ping 100
#define pong 101
#define number_of_messages 50
#define length_of_message  10

int main(int argc, char* argv[])
{
	MPI_ERR_CHECK(MPI_Init(&argc, &argv));

	int my_rank;
	MPI_ERR_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &my_rank));

	// Fill sample buffer.
	int buffer[length_of_message];
	for (int i = 0; i < length_of_message; i++)
		buffer[i] = i;

	if (my_rank == proc_A) 
	{
		double start = MPI_Wtime();

		for (int i = 0; i < number_of_messages; i++)
		{
			// TODO Send my current buffer content to proc_B (ping)

			// TODO Receive proc_B's response (pong) into my data buffer
		}

		double finish = MPI_Wtime();

		double time = finish - start;

		printf("Time for one messsage: %f seconds.\n", 
			(float)(time / (2 * number_of_messages)));
	}
	else if (my_rank == proc_B) 
	{
		for (int i = 0; i < number_of_messages; i++)
		{
			// TODO Receive proc_A's ping into my data buffer

			// Slightly modify the data
			for (int i = 0; i < length_of_message; i++)
				buffer[i]++;

			// TODO Send my current buffer content back to proc_A (pong)
		}
	} 

	// Check results correctness.
        if (my_rank == proc_A) {
		for (int i = 0; i < length_of_message; i++)
			if (buffer[i] != i + number_of_messages )
				printf("Error @ i = %d: %d != %d\n", i, buffer[i], i + 1);
        }

	MPI_ERR_CHECK(MPI_Finalize());

	return 0;
}

