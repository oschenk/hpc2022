from mpi4py import MPI
import numpy as np
import sys
import time
from pde_miniapp_py import data, linalg, operators

def readcmdline():
    """
    Read command line arguments
    """
    # report usage
    if len(sys.argv) < 4 or len(sys.argv) > 5:
        print("Usage: main nx ny nt t")
        print("  nx  number of gridpoints in x- and y-direction, respectively")
        print("  nt  number of time steps")
        print("  tf  final time")
        print("  v   (optional) turn on verbose output and write every step")
        sys.exit()

    # read nx
    nx = int(sys.argv[1])
    if nx < 1:
        sys.exit("nx must be a positive integer")
    ny = nx

    # read nt
    nt = int(sys.argv[2])
    if nt < 1:
        sys.exit("nt must be a positive integer")

    # read final time
    tf = float(sys.argv[3])
    if tf < 0:
        sys.exit("tf must be a positive real value")

    # verbose output
    verbose_output = False
    if len(sys.argv) == 5:
        verbose_output = True

    # compute time step size
    dt = tf / nt

    # compute the distance between grid points
    # assume that x/y dimension has length 1.0
    dx = 1. / (nx - 1)
    dy = dx

    # set alpha
    D = 1. # diffusion coefficient
    alpha = dx**2 / (D * dt)

    #  set beta
    R = 500. # reaction coefficient
    beta = (R * dx**2) / D
    

    return data.Discretization(nx, ny, nt, dx, dy, dt, alpha, beta), \
           verbose_output

def write_output(fname, field, timestep, time, timetable="timetable.txt",
                 hdf5=False):
    """
    Write field to output file (fname_timestep.dat) and append time
    to time table (timetable.txt)
    """
    # write field
    if not hdf5:
        field.write_mpiio(f"{fname}_{timestep:04d}.dat")
    else:
        field.write_phdf5(f"{fname}_{timestep:04d}.h5")

    # write time table (only master process!)
    if field.domain.rank == 0:
        with open(timetable, "a") as f:
            if timestep == 0:
                f.seek(0)
                f.truncate()
                f.write("# {:4s} {:10s}\n".format("step", "time"))
            f.write(f"  {timestep:4d} {time:10f}\n")

if __name__ == "__main__":
    # read cmmand line arguments
    options, verbose_output = readcmdline()

    # get COMM_WORLD communicator, size & rank
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    # initialize (sub)domain
    domain = data.Domain(comm, options)

    # set some parameters
    max_cg_iters = 200
    max_newton_iters = 50
    tolerance = 1.e-6
    skip = 10 # write output every skip time step (or just last if skip < 1)

    if rank == 0:
        print(80*"=")
        print('{:^80}'.format("Welcome to mini-stencil!"))
        print("version   :: MPI Python")
        print(f"mesh      :: {options.nx:4d} x {options.ny:4d},",
              f"dx = {options.dx:f}, dy = {options.dy:f}")
        print("time      :: {nt:4d} time steps from 0 to {time:f}"
             .format(nt=options.nt, time=options.nt*options.dt))
        print(f"iteration :: CG {max_cg_iters:d},",
                            f"Newton {max_newton_iters:d},",
                            f"tolerance {tolerance:e}")
        print(80*"=")

    # print some domain decomposition specifics
    domain.print()

    # allocate global fields
    s_old  = data.Field(domain)
    s_new  = data.Field(domain)
    f      = data.Field(domain)
    deltas = data.Field(domain)
    v      = data.Field(domain)
    fv     = data.Field(domain)
    Fx     = data.Field(domain)

    # set the initial condition (boundary conditions are already set to 0)
    xc = 1. / 4.
    yc = 1. / 4.
    radius = 1. / 8.
    s_new.inner[...] = 0.1
    s_new.bdryE[...] = s_new.bdryN[...] = s_new.bdryS[...] = s_new.bdryW[...] = 0.1
    s_old.bdryE[...] = s_old.bdryN[...] = s_old.bdryS[...] = s_old.bdryW[...] = 0.1
    v.bdryE[...] = v.bdryN[...] = v.bdryS[...] = v.bdryW[...] = 0.1

    x_global  = (  s_new.domain.global_startx
                 + np.arange(s_new.domain.local_nx))*options.dx
    y_global  = (  s_new.domain.global_starty
                 + np.arange(s_new.domain.local_ny))*options.dy
    X_global, Y_global = np.meshgrid(x_global, y_global,
                                     indexing='ij', sparse=True)
    R2_global = (X_global - xc)**2 + (Y_global - yc)**2
    s_new.inner[R2_global < radius**2] = 0.2

    # write initial conditions
    write_output("simulation", s_new, 0, 0.)

    # set some counters
    iters_newton = 0
    iters_cg     = 0

    # initialize CG solver
    cg = linalg.hpc_cg(s_new.domain)

    if domain.rank == 0 and verbose_output:
        print("=== Begin of simulation ===")
        print(80*"-")

    # start timer
    timespent = - time.perf_counter()

    # main time loop
    t        = 0. # initial time
    timestep = 0 # initial step
    while timestep < options.nt:
        # print some info
        if domain.rank == 0 and verbose_output:
            print(f"timestep / time: {timestep:5d} / {t:10f}")

        # set s_new and s_old to be the solution
        s_old.inner[...] = s_new.inner[...]

        converged = False
        for it in range(max_newton_iters):
            # compute residual
            operators.diffusion(options, s_old, s_new, f)
            residual = linalg.hpc_norm2(f)

            # check for convergence
            if residual < tolerance:
                converged = True
                break

            def A(x, Ax):
                eps = 1.e-4
                v.inner[...] = s_new.inner[...] + eps*x.inner[...]
                # f(v = s_new + eps*x)
                operators.diffusion(options, s_old, v, Fx)
                # Df(s_new)*x = A*x ~= (fv - f)/eps
                Ax.inner[...] = (Fx.inner[...] - f.inner[...])/eps

            # solve linear system to get -deltax
            cg_converged, it_cg = cg.solve(A, f, deltas, tolerance,
                                           max_cg_iters)
            iters_cg += it_cg

            # check that the CG solver converged
            if not cg_converged:
                break

            # update solution
            s_new.inner[...] -= deltas.inner[...]

        # check result of Newton iterations
        if converged:
            timestep += 1
            t        += options.dt
            iters_newton += it + 1
            # write output
            if (skip > 0 and timestep % skip == 0) or timestep == options.nt:
                write_output("simulation", s_new, timestep, t)
            if domain.rank == 0:
                if verbose_output:
                    print(f"success, new time: {t:10f}")
                    print(f"required {it:5d} Newton iterations",
                          f"for residual {residual:10e}")
                    print(80*"-")
        else:
            if domain.rank == 0:
                print(f"step {timestep:d} ERROR:",
                       "Newton iterations failed to converge,",
                      f"residual = {residual:f}")
            break


    # stop timer
    timespent += time.perf_counter()

    # print table sumarizing results
    if domain.rank == 0:
        print(f"simulation took {timespent:f} seconds")
        print(f"{iters_cg:d} conjugate gradient iterations:",
               "at rate of {:f} iters/second".format(iters_cg/timespent))
        print(f"{iters_newton:d} newton iterations")
        print(80*"-")

    if domain.rank == 0:
        print("=== End of simulation. Goodbye! ===")

