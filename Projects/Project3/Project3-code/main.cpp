// ******************************************
// implicit time stepping implementation of 2D diffusion problem
// Ben Cumming, CSCS
// *****************************************

// A small benchmark app that solves the 2D fisher equation using second-order
// finite differences.

// Syntax: ./main nx nt t

#include <algorithm>
#include <fstream>
#include <iostream>

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <omp.h>
#include <stdio.h>

#include "data.h"
#include "linalg.h"
#include "operators.h"
#include "stats.h"

using namespace data;
using namespace linalg;
using namespace operators;
using namespace stats;

// ==============================================================================

// read command line arguments
static void readcmdline(Discretization& options, int argc, char* argv[])
{
    if (argc<4 || argc>5)
    {
        printf("Usage: main nx nt t verbose\n");
        printf("  nx        number of gridpoints in x-direction and y-direction, respectively\n");
        printf("  nt        number of timesteps\n");
        printf("  t         total time\n");
        printf("  verbose   (optional) verbose output\n");
        exit(1);
    }

    // read nx
    options.nx = atoi(argv[1]);
    if (options.nx < 1)
    {
        fprintf(stderr, "nx must be positive integer\n");
        exit(-1);
    }

    // read nt
    options.nt = atoi(argv[2]);
    if (options.nt < 1)
    {
        fprintf(stderr, "nt must be positive integer\n");
        exit(-1);
    }

    // read total time
    double t = atof(argv[3]);
    if (t < 0)
    {
        fprintf(stderr, "t must be positive real value\n");
        exit(-1);
    }

    // set verbosity if requested
    verbose_output = false;
    if (argc==5) {
        verbose_output = true;
    }

    // store the parameters
    options.N = options.nx * options.nx;

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
    int nx = options.nx;
    int N  = options.N;
    int nt = options.nt;

    // set iteration parameters
    int max_cg_iters     = 200;
    int max_newton_iters = 50;
    double tolerance     = 1.e-6;

    std::cout << "========================================================================" << std::endl;
    std::cout << "                      Welcome to mini-stencil!" << std::endl;
    std::cout << "version   :: Serial C++" << std::endl;
    std::cout << "mesh      :: " << options.nx << " * " << options.nx << " dx = " << options.dx << std::endl;
    std::cout << "time      :: " << nt << " time steps from 0 .. " << options.nt*options.dt << std::endl;;
    std::cout << "iteration :: " << "CG "          << max_cg_iters
                                 << ", Newton "    << max_newton_iters
                                 << ", tolerance " << tolerance << std::endl;;
    std::cout << "========================================================================" << std::endl;

    // allocate global fields
    // allocate global fields
    y_new.init(nx,nx);
    y_old.init(nx,nx);
    bndN.init(nx,1);
    bndS.init(nx,1);
    bndE.init(nx,1);
    bndW.init(nx,1);

    Field b(nx,nx);
    Field deltax(nx,nx);

    // set dirichlet boundary conditions to 0 all around
    double const_bdy = 0.1;
    hpc_fill(bndN, const_bdy, nx);
    hpc_fill(bndS, const_bdy, nx);
    hpc_fill(bndE, const_bdy, nx);
    hpc_fill(bndW, const_bdy, nx);

    // set the initial condition
    // a circle of concentration 0.1 centred at (xdim/4, ydim/4) with radius
    // no larger than 1/8 of both xdim and ydim
    double const_fill = 0.1;
    double inner_circle = 0.2;
    hpc_fill(y_new, const_fill, nx*nx);
    double xc = 1.0 / 4.0;
    double yc = (nx - 1) * options.dx / 4;
    double radius = fmin(xc, yc) / 2.0;
    for (int j = 0; j < nx; j++)
    {
        double y = (j - 1) * options.dx;
        for (int i = 0; i < nx; i++)
        {
            double x = (i - 1) * options.dx;
            if ((x - xc) * (x - xc) + (y - yc) * (y - yc) < radius * radius)
                y_new[i+nx*j] = inner_circle;
        }
    }

    iters_cg = 0;
    iters_newton = 0;

    // start timer
    double timespent = -omp_get_wtime();

    // main timeloop
    for (int timestep = 1; timestep <= nt; timestep++)
    {
        // set y_new and y_old to be the solution
        hpc_copy(y_old, y_new, N);

        double residual;
        bool converged = false;
        int it;
        for (it=0; it<max_newton_iters; it++)
        {
            // compute residual : requires both y_new and y_old
            diffusion(y_new, b);
            residual = hpc_norm2(b, N);

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
            hpc_axpy(y_new, -1.0, deltax, N);
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
    timespent += omp_get_wtime();

    ////////////////////////////////////////////////////////////////////
    // write final solution to BOV file for visualization
    ////////////////////////////////////////////////////////////////////

    // binary data
    {
        FILE* output = fopen("output.bin", "w");
        fwrite(y_new.data(), sizeof(double), nx * nx, output);
        fclose(output);
    }

    std::ofstream fid("output.bov");
    fid << "TIME: 0.0" << std::endl;
    fid << "DATA_FILE: output.bin" << std::endl;
    fid << "DATA_SIZE: " << options.nx << " " << options.nx << " 1" << std::endl;;
    fid << "DATA_FORMAT: DOUBLE" << std::endl;
    fid << "VARIABLE: phi" << std::endl;
    fid << "DATA_ENDIAN: LITTLE" << std::endl;
    fid << "CENTERING: nodal" << std::endl;
    fid << "BRICK_SIZE: 1.0 " << (options.nx-1)*options.dx << " 1.0" << std::endl;

    // print table sumarizing results
    std::cout << "--------------------------------------------------------------------------------"
              << std::endl;
    std::cout << "simulation took " << timespent << " seconds" << std::endl;
    std::cout << int(iters_cg) << " conjugate gradient iterations, at rate of "
              << float(iters_cg)/timespent << " iters/second" << std::endl;
    std::cout << iters_newton << " newton iterations" << std::endl;
    std::cout << "--------------------------------------------------------------------------------"
              << std::endl;

    std::cout << "Goodbye!" << std::endl;

    return 0;
}

