"""Collection of visualization tools"""
import matplotlib.pyplot as plt
import numpy as np
import struct

def read_time(timetable):
    """Read time table"""
    return np.reshape(np.loadtxt(timetable), (-1, 2))

def read_field(fname):
    """Read data"""
    with open(fname, "rb") as f:
        # read header
        nx = struct.unpack("i", f.read(4))[0]
        ny = struct.unpack("i", f.read(4))[0]

        # read data
        tmp = np.fromfile(f, dtype=np.double)

    return tmp.reshape((nx, ny))

def draw(model, step, timetable="timetable.txt"):
    """Draw model output at step"""

    # get field
    A = read_field(f"{model}_{step:04d}.dat")
    nx = A.shape[0]
    ny = A.shape[1]

    # get time
    if timetable != None:
        tt = read_time(timetable)
        time = np.asscalar(tt[np.where(tt[:,0].astype(int) == step), 1])
    else:
        time = None

    # create grid
    x = np.linspace(0., 1., nx)
    y = np.linspace(0., 1., ny)
    [X, Y] = np.meshgrid(x, y)

    # draw
    plt.clf()
    N = 12 # number of contours
    C1 = plt.contourf(X, Y, A.T, N, alpha=0.75, cmap="jet")
    C2 = plt.contour (X, Y, A.T, N, colors='black', linewidths=0.1)
    

    # plot cosmetics
    plt.gcf().colorbar(C1)
    plt.clabel(C2, inline=1)
    plt.axis("scaled")
    plt.xlabel("x")
    plt.xlabel("y")
    if time != None:
        title = f"step / time: {step:4d} / {time:10f}"
    else:
        title = f"step : {step:4d}"
    plt.title(title)

