from mpi4py import MPI

# get comm, size & rank
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

# hello
print(f"Hello world from rank {rank} out of {size} processes")

