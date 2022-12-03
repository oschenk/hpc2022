"""Collection of data structures"""
import h5py
import math
from mpi4py import MPI
import numpy as np
import struct
import time

class Discretization:
    """Simple class to hold discretization specifics"""
    def __init__(self, nx, ny, nt, dx, dy, dt, alpha, beta):
        self._nx = nx
        self._ny = ny
        self._nt = nt
        self._dx = dx
        self._dy = dy
        self._dt = dt
        self._alpha = alpha
        self._beta  = beta

    @property
    def nx(self):
        """Number of gridpoints in x-direction"""
        return self._nx

    @property
    def ny(self):
        """Number of gridpoints in y-direction"""
        return self._ny

    @property
    def nt(self):
        """Number of time steps"""
        return self._nt

    @property
    def dt(self):
        """Time step"""
        return self._dt

    @property
    def dx(self):
        """Grid spacing in x-direction"""
        return self._dx

    @property
    def dy(self):
        """Grid spacing in y-direction"""
        return self._dy

    @property
    def alpha(self):
        """alpha = dx**2 / (D * dt)"""
        return self._alpha

    @property
    def beta(self):
        """beta = (R * dx**2)/D"""
        return self._beta

class Domain:
    """Simple class to hold two-dimensional (sub)domain specifics"""
    def __init__(self, comm, discretization):
        # set comm
        self._comm = comm

        # get size & rank
        self._size = comm.Get_size()
        self._rank = comm.Get_rank()

        # compute a distribution of processes per coordinate direction
        self._dims = [0, 0]
        self._dims = MPI.Compute_dims(self._size, self._dims)

        # create two-dimensional non-periodic Cartesian topology
        self._periods = [False, False]
        self._comm_cart = self._comm.Create_cart(self._dims,
                                                 periods=self._periods)

        # get rank's coordinates in the topology
        self._coords = self._comm_cart.Get_coords(self._rank)

        # get rank's south/north/west/east neighbors
        self._neigh_west,  self._neigh_east  = self._comm_cart.Shift(0, 1)
        self._neigh_south, self._neigh_north = self._comm_cart.Shift(1, 1)

        # global domain discretization size
        self._global_nx = discretization.nx
        self._global_ny = discretization.ny

        # global start/end index in (sub)domain discretization
        chunkx = discretization.nx // self._dims[0]
        self._global_startx = chunkx * self._coords[0]
        if self._coords[0] == self._dims[0] - 1:
            self._global_endx = discretization.nx
        else:
            self._global_endx = self._global_startx + chunkx
        chunky = discretization.ny // self._dims[1]
        self._global_starty = chunky * self._coords[1]
        if self._coords[1] == self._dims[1] - 1:
            self._global_endy = discretization.ny
        else:
            self._global_endy = self._global_starty + chunky

        # local (sub)domain discretization size
        self._local_nx = self._global_endx - self._global_startx
        self._local_ny = self._global_endy - self._global_starty

    @property
    def comm(self):
        """MPI communicator"""
        return self._comm

    @property
    def size(self):
        """MPI size"""
        return self._size

    @property
    def rank(self):
        """MPI rank"""
        return self._rank

    @property
    def dims(self):
        """"Number of processes in x and y dimension"""
        return self._dims

    @property
    def comm_cart(self):
        """MPI Cartesian topology communicator"""
        return self._comm_cart

    @property
    def coords(self):
        """Rank's Cartesian coordinates in Cartesian topology"""
        return self._coords

    @property
    def neighbour_north(self):
        """Rank of neighbouring process in North direction"""
        return self._neigh_north

    @property
    def neighbour_east(self):
        """Rank of neighbouring process in East direction"""
        return self._neigh_east

    @property
    def neighbour_south(self):
        """Rank of neighbouring process in South direction"""
        return self._neigh_south

    @property
    def neighbour_west(self):
        """Rank of neighbouring process in West direction"""
        return self._neigh_west

    @property
    def global_startx(self):
        """Return global start index of (sub)domain in x-direction"""
        return self._global_startx

    @property
    def global_starty(self):
        """Return global start index of (sub)domain in y-direction"""
        return self._global_starty

    @property
    def global_endx(self):
        """Return global end index of (sub)domain in x-direction"""
        return self._global_endx

    @property
    def global_endy(self):
        """Return global end index of (sub)domain in y-direction"""
        return self._global_endy

    @property
    def global_nx(self):
        """Return global (sub)domain size in x-direction"""
        return self._global_nx

    @property
    def global_ny(self):
        """Return global (sub)domain size in y-direction"""
        return self._global_ny

    @property
    def local_nx(self):
        """Return local (sub)domain size in x-direction"""
        return self._local_nx

    @property
    def local_ny(self):
        """Return local (sub)domain size in y-direction"""
        return self._local_ny

    def print(self):
        """Print domain decomposition specifics"""
        if self._rank == 0:
            print("{:d} processes decomposed into a".format(self._size),
                  "{:d} x {:d} Cartesian topology".format(self._dims[0],
                                                          self._dims[1]),
                  flush=True)
        self._comm.Barrier()
        for irank in range(self._size):
            self._comm.Barrier()
            if self._rank == irank:
                print(MPI.Get_processor_name(),
                      "rank {:3d} / {:3d} :".format(self._rank, self._size),
                      "({:3d},{:3d})".format(self._coords[0], self._coords[1]),
                      "neigh N:S {:3d} : {:3d}".format(self._neigh_north,
                                                       self._neigh_south),
                      "neigh E:W {:3d} : {:3d}".format(self._neigh_east,
                                                       self._neigh_west),
                      "local size {:3d} x {:3d}".format(self._local_nx,
                                                        self._local_ny),
                      flush=True)

        # add artificial pause so that output doesn't pollute later output
        time.sleep(0.1)

class Field:
    """Simple class to hold a distributed field"""
    def __init__(self, domain):
        # set Domain
        self._domain = domain

        # inner field data
        self._inner = np.zeros((domain.local_nx, domain.local_ny), dtype='d')

        # north/south/east/west boundary field data
        self._bdryN  = np.zeros(domain.local_nx, dtype='d')
        self._bdryE  = np.zeros(domain.local_ny, dtype='d')
        self._bdryS  = np.zeros(domain.local_nx, dtype='d')
        self._bdryW  = np.zeros(domain.local_ny, dtype='d')

        # north/south/east/west boundary buffer data for MPI
        self._buffN  = np.zeros(domain.local_nx, dtype='d')
        self._buffE  = np.zeros(domain.local_ny, dtype='d')
        self._buffS  = np.zeros(domain.local_nx, dtype='d')
        self._buffW  = np.zeros(domain.local_ny, dtype='d')

    @property
    def domain(self):
        """(Sub)domain of distributed field"""
        return self._domain

    @property
    def inner(self):
        """Inner part of distributed field"""
        return self._inner

    @property
    def bdryN(self):
        """North boundary data of distributed field"""
        return self._bdryN

    @property
    def bdryE(self):
        """East boundary data of distributed field"""
        return self._bdryE

    @property
    def bdryS(self):
        """South boundary data of distributed field"""
        return self._bdryS

    @property
    def bdryW(self):
        """West boundary data of distributed field"""
        return self._bdryW

    def exchange_startall(self):
        """Start exchanging boundary field data"""
        domain = self._domain # copy for convenience
        # ... implement ...

    def exchange_waitall(self):
        """Wait until exchanging boundary field data is complete"""
        # ... implement ...

    def write_mpiio(self, fname):
        """Write field to file fname with MPI-IO"""

        # open file
        amode = MPI.MODE_CREATE | MPI.MODE_WRONLY
        fh = MPI.File.Open(self._domain.comm_cart, fname, amode)

        # get numpy data type and convert to MPI data type
        # (note that we are using a numpy "private" attribute)
        mpi_t = MPI._typedict[self._inner.dtype.char]

        # create file type
        data_sizes    = [self._domain.global_nx,
                         self._domain.global_ny]
        data_subsizes = [self._domain.local_nx,
                         self._domain.local_ny ]
        data_starts   = [self._domain.global_startx,
                         self._domain.global_starty]
        file_t = mpi_t.Create_subarray(data_sizes, data_subsizes, data_starts)
        file_t.Commit()

        # set file view
        fh.Set_view(8, filetype=file_t)

        # write data
        fh.Write_all(self._inner)

        # free file type
        file_t.Free()

        # close file
        fh.Close()

        # write header
        if self._domain.rank == 0:
            with open(fname, "r+b") as f:
                f.write(struct.pack('i', data_sizes[0]))
                f.write(struct.pack('i', data_sizes[1]))

    def write_phdf5(self, fname, field_name="field"):
        """Write field to file fname with parallel HDF5"""

        # open file
        f = h5py.File(fname, "w", driver="mpio", comm=self._domain.comm_cart)

        # create data set
        data_sizes  = (self._domain.global_nx,
                       self._domain.global_ny)
        dset = f.create_dataset(field_name, data_sizes,
                                dtype=self._inner.dtype.char)

        # write data
        starts = [self._domain.global_startx,
                  self._domain.global_starty]
        ends   = [self._domain.global_endx,
                  self._domain.global_endy]
        dset[starts[0]:ends[0], starts[1]:ends[1]] = self._inner

        # close file
        f.close()

if __name__ == '__main__':
    pass

