// ******************************************
// implicit time stepping implementation of 2D diffusion problem
// Ben Cumming, CSCS
// *****************************************

// A small benchmark app that solves the 2D fisher equation using second-order
// finite differences.

// Syntax: ./main nx nt t

#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <mpi.h>

#include "data.h"
#include "linalg.h"
#include "operators.h"
#include "stats.h"

using namespace data;
using namespace linalg;
using namespace operators;
using namespace stats;

// ==============================================================================
void write_binary(std::string fname, Field &u, SubDomain &domain, Discretization &options)
{
    MPI_Offset disp = 0;
    MPI_File filehandle;
    MPI_Datatype filetype;

    int result =
        MPI_File_open(
            domain.comm_cart,
            fname.c_str(),
            MPI_MODE_CREATE | MPI_MODE_WRONLY,
            MPI_INFO_NULL,
            &filehandle
        );
    assert(result==MPI_SUCCESS);

    int ustart[]  = {domain.startx-1, domain.starty-1};
    int ucount[]  = {domain.nx, domain.ny};
    int dimuids[] = {options.nx, options.ny};

    result = MPI_Type_create_subarray(2, dimuids, ucount, ustart, MPI_ORDER_FORTRAN, MPI_DOUBLE, &filetype);
    assert(result==MPI_SUCCESS);

    result = MPI_Type_commit(&filetype);
    assert(result==MPI_SUCCESS);

    result = MPI_File_set_view(filehandle, disp, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
    assert(result==MPI_SUCCESS);

    result = MPI_File_write_all(filehandle, u.data(), domain.N, MPI_DOUBLE, MPI_STATUS_IGNORE);
    assert(result==MPI_SUCCESS);

    result = MPI_Type_free(&filetype);
    assert(result==MPI_SUCCESS);

    result = MPI_File_close(&filehandle);
    assert(result==MPI_SUCCESS);
}

// read command line arguments
static void readcmdline(Discretization& options, int argc, char* argv[])
{
    if (argc<4 || argc>5)
    {
        std::cerr << "Usage: main nx nt t\n";
        std::cerr << "  nx  number of gridpoints in x-direction and y-direction, respectively\n";
        std::cerr << "  nt  number of timesteps\n";
        std::cerr << "  t   total time\n";
        std::cerr << "  v   [optional] turn on verbose output\n";
        exit(1);
    }

    // read nx
    options.nx = atoi(argv[1]);
    if (options.nx < 1) {
        std::cerr << "nx must be positive integer\n";
        exit(-1);
    }
    options.ny = options.nx;

    // read nt
    options.nt = atoi(argv[2]);
    if (options.nt < 1)
    {
        std::cerr << "nt must be positive integer\n";
        exit(-1);
    }

    // read total time
    double t = atof(argv[3]);
    if (t < 0)
    {
        std::cerr << "t must be positive real value\n";
        exit(-1);
    }

    // set verbosity if requested
    verbose_output = false;
    if (argc==5) {
        verbose_output = (domain.rank==0);
    }

    // compute timestep size
    options.dt = t / options.nt;

    // compute the distance between grid points
    // assume that x dimension has length 1.0
    options.dx = 1. / (options.nx - 1);

    // set alpha, assume diffusion coefficient D is 1
    options.alpha = (options.dx * options.dx) / (1. * options.dt);

    // set beta, assume diffusion coefficient D=1, reaction coefficient R=1000
    double R = 500.;
    options.beta = (R * options.dx * options.dx)/1.;
}

// ==============================================================================

int main(int argc, char* argv[])
{
    // read command line arguments
    readcmdline(options, argc, argv);

    // initialize MPI
    int mpi_rank, mpi_size, threadLevelProvided;
    // TODO initialize
    // use "MPI_Comm_size", "MPI_Comm_rank" and "MPI_Init_thread"

    // initialize subdomain
    domain.init(mpi_rank, mpi_size, options);
    domain.print();

    int nx = domain.nx;
    int ny = domain.ny;
    int N  = domain.N;
    int nt  = options.nt;

    // set iteration parameters
    int max_cg_iters     = 200;
    int max_newton_iters = 50;
    double tolerance     = 1.e-6;

    if( domain.rank == 0 ) {
        std::cout << "========================================================================" << std::endl;
        std::cout << "                      Welcome to mini-stencil!" << std::endl;
        std::cout << "version   :: MPI : " << domain.size << " MPI ranks" << std::endl;
        std::cout << "mesh      :: " << options.nx << " * " << options.ny << " dx = " << options.dx << std::endl;
        std::cout << "time      :: " << nt << " time steps from 0 .. " << options.nt*options.dt << std::endl;;
        std::cout << "iteration :: " << "CG "          << max_cg_iters
                                     << ", Newton "    << max_newton_iters
                                     << ", tolerance " << tolerance << std::endl;;
        std::cout << "========================================================================" << std::endl;
    }

    // allocate global fields
    y_new.init(nx,ny);
    y_old.init(nx,ny);
    bndN.init(nx,1);
    bndS.init(nx,1);
    bndE.init(ny,1);
    bndW.init(ny,1);
    buffN.init(nx,1);
    buffS.init(nx,1);
    buffE.init(ny,1);
    buffW.init(ny,1);

    Field b(nx,ny);
    Field deltax(nx,ny);

    // set dirichlet boundary conditions to 0.1 all around
    double bdy_value = 0.1;
    hpc_fill(bndN, bdy_value);
    hpc_fill(bndS, bdy_value);
    hpc_fill(bndE, bdy_value);
    hpc_fill(bndW, bdy_value);

    // set the initial condition
    // a circle of concentration 0.1 centred at (xdim/4, ydim/4) with radius
    // no larger than 1/8 of both xdim and ydim
    double const_fill = 0.1;
    double inner_circle = 0.2;
    hpc_fill(y_new, const_fill);
    double xc = 1.0 / 4.0;
    double yc = (options.ny - 1) * options.dx / 4;
    double radius = std::min(xc, yc) / 2.0;
    for (int j = domain.starty-1; j < domain.endy; j++)
    {
        double y = (j - 1) * options.dx;
        for (int i = domain.startx-1; i < domain.endx; i++)
        {
            double x = (i - 1) * options.dx;
            if ((x - xc) * (x - xc) + (y - yc) * (y - yc) < radius * radius)
                y_new(i-domain.startx+1, j-domain.starty+1) = inner_circle;
        }
    }

    iters_cg = 0;
    iters_newton = 0;

    // start timer
    double timespent = -MPI_Wtime();

    // main timeloop
    for (int timestep = 1; timestep <= nt; timestep++)
    {
        // set y_new and y_old to be the solution
        hpc_copy(y_old, y_new);

        double residual;
        bool converged = false;
        int it;
        for (it=0; it<max_newton_iters; it++)
        {
            // compute residual : requires both y_new and y_old
            diffusion(y_new, b);
            residual = hpc_norm2(b);

            // check for convergence
            if (residual < tolerance)
            {
                converged = true;
                break;
            }

            // solve linear system to get -deltax
            bool cg_converged = false;
            hpc_cg(deltax, y_new, b, max_cg_iters, tolerance, cg_converged);

            // check that the CG solver converged
            if (!cg_converged) break;

            // update solution
            hpc_axpy(y_new, -1.0, deltax);
        }
        iters_newton += it+1;

        // output some statistics
        if (converged && verbose_output) {
            std::cout << "step " << timestep
                      << " required " << it
                      << " iterations for residual " << residual
                      << std::endl;
        }
        if (!converged) {
            std::cerr << "step " << timestep
                      << " ERROR : nonlinear iterations failed to converge" << std::endl;;
            break;
        }
    }

    // get times
    timespent += MPI_Wtime();

    ////////////////////////////////////////////////////////////////////
    // write final solution to BOV file for visualization
    ////////////////////////////////////////////////////////////////////

    // binary data
    write_binary("output.bin", y_old, domain, options);

    // metadata
    if( domain.rank==0 ) {
        std::ofstream fid("output.bov");
        fid << "TIME: 0.0" << std::endl;
        fid << "DATA_FILE: output.bin" << std::endl;
        fid << "DATA_SIZE: " << options.nx << " " << options.ny << " 1" << std::endl;;
        fid << "DATA_FORMAT: DOUBLE" << std::endl;
        fid << "VARIABLE: phi" << std::endl;
        fid << "DATA_ENDIAN: LITTLE" << std::endl;
        fid << "CENTERING: nodal" << std::endl;
        fid << "BRICK_SIZE: 1.0 " << (options.ny-1)*options.dx << " 1.0" << std::endl;
    }

    // print table sumarizing results
    if(domain.rank == 0) {
        std::cout << "--------------------------------------------------------------------------------"
                  << std::endl;
        std::cout << "simulation took " << timespent << " seconds" << std::endl;
        std::cout << int(iters_cg) << " conjugate gradient iterations, at rate of "
                  << float(iters_cg)/timespent << " iters/second" << std::endl;
        std::cout << iters_newton << " newton iterations" << std::endl;
        std::cout << "--------------------------------------------------------------------------------"
                  << std::endl;
    }

    if(domain.rank==0)
        std::cout << "Goodbye!" << std::endl;

    // TODO finalize it using "MPI_Finalize" and "MPI_Comm_free"

    return 0;
}

