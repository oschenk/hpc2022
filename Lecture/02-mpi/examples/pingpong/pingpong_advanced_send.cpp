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
 * Purpose: Benchmarking MPI_Send and MPI_Recv.                 *
 *                                                              *
 * Contents: C-Source                                           *
 *                                                              *
 ****************************************************************/

#include <cstdio>
#include <cstdlib>
#include <mpi.h>
#include <vector>

#define MPI_ERR_CHECK(call)                          \
    do { int err = call; if (err != MPI_SUCCESS) {   \
        char errstr[MPI_MAX_ERROR_STRING];           \
        int szerrstr;                                \
        MPI_Error_string(err, errstr, &szerrstr);    \
        fprintf(stderr, "MPI error at %s:%i : %s\n", \
            __FILE__, __LINE__, errstr);             \
        MPI_Abort(MPI_COMM_WORLD, 1);                \
    }} while (0);

#define process_A 0
#define process_B 1
#define ping 100
#define pong 101

#define number_of_messages 50
#define start_length 8
#define length_factor 64
#define max_length 2097152     /* 2 Mega */
#define number_package_sizes 4

using namespace std;

int main(int argc, char* argv[])
{
	MPI_ERR_CHECK(MPI_Init(&argc, &argv));

	int my_rank;
	MPI_ERR_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &my_rank));

	int length_of_message;
	vector<float> vbuffer;
	vbuffer.resize(max_length);
	float* buffer = &vbuffer[0];
	if (my_rank == process_A)
	{
		length_of_message = start_length;

		printf("message size\ttransfer time\t\tbandwidth\n");

		for (int i = 0; i < number_package_sizes; i++)
		{ 
			double start = MPI_Wtime();

			for (int j = 0; j < number_of_messages; j++)
			{
				MPI_ERR_CHECK(MPI_Send(buffer, length_of_message, MPI_FLOAT, process_B, 
					ping, MPI_COMM_WORLD));

				MPI_Status status;
				MPI_ERR_CHECK(MPI_Recv(buffer, length_of_message, MPI_FLOAT, process_B, 
					pong, MPI_COMM_WORLD, &status));
			}

			double finish = MPI_Wtime();

			double time = finish - start;

			double transfer_time = time / (2 * number_of_messages);

			printf("%zu bytes\t\t%f sec\t\t%f MB/s\n", 
				length_of_message * sizeof(float), transfer_time,
				1.0e-6 * length_of_message * sizeof(float) / transfer_time);

			length_of_message *= length_factor;
		}
	}
	else if (my_rank == process_B) 
	{
		length_of_message = start_length;

		for (int i = 0; i < number_package_sizes; i++)
		{    
			for (int j = 0; j < number_of_messages; j++)
			{
				MPI_Status status;
				MPI_ERR_CHECK(MPI_Recv(buffer, length_of_message, MPI_FLOAT, process_A, 
					ping, MPI_COMM_WORLD, &status));

				MPI_ERR_CHECK(MPI_Send(buffer, length_of_message, MPI_FLOAT, process_A, 
					pong, MPI_COMM_WORLD));
			}  
			length_of_message *= length_factor;
		}
	}

	MPI_ERR_CHECK(MPI_Finalize());

	return 0;
}

