#include <stdio.h>
#include <mpi.h>

int main (int argc,char** argv) {
  int rank, size, dest_proc, from_proc;
  int message[1];
  int buffer[1];
  int tag = 0;
  MPI_Status stat;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  *message = rank;

  from_proc = rank -1;
  dest_proc = rank +1;
  if(rank == 0) {
    printf("Hello world!\n");
    from_proc = size-1;
  }
  if(rank == size-1){
    dest_proc = 0;
  }
  // begin snippet
  if(size > 1) {
    if(rank == 0) {
      printf("This is process %d out of %d.\n", rank, size);
      MPI_Send(message, 1, MPI_INT, dest_proc, tag, MPI_COMM_WORLD);
      MPI_Recv(buffer, 1,MPI_INT, from_proc, tag, MPI_COMM_WORLD, &stat);
    }
    else{
      MPI_Recv(buffer, 1,MPI_INT, from_proc, tag, MPI_COMM_WORLD, &stat);
      printf("This is process %d out of %d.\n", rank, size);
      MPI_Send(message, 1, MPI_INT, dest_proc, tag, MPI_COMM_WORLD);
    }
       
  }
  else{
    printf("This is process %d out of %d.\n", rank, size);
  }
  // end snippet
  MPI_Finalize();
  return 0;
}
