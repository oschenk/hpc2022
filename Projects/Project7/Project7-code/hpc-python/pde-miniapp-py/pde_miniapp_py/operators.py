"""
Simple grid operators
"""

def diffusion(discretization, s_old, s_new, f):
    """Diffusion operator"""
    alpha = discretization.alpha
    beta  = discretization.beta

    # start exchanging
    s_new.exchange_startall()

    # work on interior points
    f.inner[1:-1, 1:-1] = (
                           # new central point
                           - (4. + alpha) * s_new.inner[1:-1, 1:-1]
                           # new west    point
                           +                s_new.inner[ :-2, 1:-1]
                           # new east    point
                           +                s_new.inner[2:  , 1:-1]
                           # new north   point
                           +                s_new.inner[1:-1, 2:  ]
                           # new south   point
                           +                s_new.inner[1:-1,  :-2]
                           # old central point
                           +       alpha  * s_old.inner[1:-1, 1:-1]
                           # new central reaction
                           + beta *        s_new.inner[1:-1, 1:-1]
                                  * (1.0 - s_new.inner[1:-1, 1:-1])
                          )
    # wait on echanges to complete
    s_new.exchange_waitall()

    # work on west and east boundary
    for i in [0, -1]:
        if   i ==  0: # west boundary
            s_new_west = s_new.bdryW[1:-1]
            s_new_east = s_new.inner[i+1, 1:-1]
        elif i == -1: # east boundary
            s_new_west = s_new.inner[i-1, 1:-1]
            s_new_east = s_new.bdryE[1:-1]
        f.inner[i, 1:-1] = (
                               # new central point
                               - (4. + alpha) * s_new.inner[i, 1:-1]
                               # new west    point
                               +                s_new_west
                               # new east    point
                               +                s_new_east
                               # new north   point
                               +                s_new.inner[i, 2:  ]
                               # new south   point
                               +                s_new.inner[i,  :-2]
                               # old central point
                               +       alpha  * s_old.inner[i, 1:-1]
                               # new central reaction
                               + beta *        s_new.inner[i, 1:-1]
                                      * (1.0 - s_new.inner[i, 1:-1])
                           )


    # work on south and north boundary
    for j in [0, -1]:
        if   j ==  0: # south boundary
            s_new_south = s_new.bdryS[1:-1]
            s_new_north = s_new.inner[1:-1, j+1]
        elif j == -1: # north boundary
            s_new_south = s_new.inner[1:-1, j-1]
            s_new_north = s_new.bdryN[1:-1]
        f.inner[1:-1, j] = (
                            # new central point
                            - (4. + alpha) * s_new.inner[1:-1, j]
                            # new west    point
                            +                s_new.inner[ :-2, j]
                            # new east    point
                            +                s_new.inner[2:  , j]
                            # new north   point
                            +                s_new_north
                            # new south   point
                            +                s_new_south
                            # old central point
                            +       alpha  * s_old.inner[1:-1, j]
                            # new central reaction
                            + beta *        s_new.inner[1:-1, j]
                                   * (1.0 - s_new.inner[1:-1, j])
                           )

    # work on NE/SE/SW/NW corners
    for i, j in [[-1, -1], [-1, 0], [0, 0], [0, -1]]:
        if   i ==  0: # west boundary
            s_new_west  = s_new.bdryW[j]
            s_new_east  = s_new.inner[i+1, j]
        elif i == -1: # east boundary
            s_new_west  = s_new.inner[i-1, j]
            s_new_east  = s_new.bdryE[j]
        if   j ==  0: # south boundary
            s_new_south = s_new.bdryS[i]
            s_new_north = s_new.inner[i, j+1]
        elif j == -1: # north boundary
            s_new_south = s_new.inner[i, j-1]
            s_new_north = s_new.bdryN[i]
        f.inner[i,j] = (
                        # new central point
                        - (4. + alpha) * s_new.inner[i, j]
                        # new west    point
                        + s_new_west
                        # new east    point
                        + s_new_east
                        # new north   point
                        + s_new_north
                        # new south   point
                        + s_new_south
                        # old central point
                        +       alpha  * s_old.inner[i, j]
                        # new central reaction
                        + beta *        s_new.inner[i, j]
                               * (1.0 - s_new.inner[i, j])
                       )

