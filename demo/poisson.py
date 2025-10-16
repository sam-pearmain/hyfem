import os 
os.environ["OMP_NUM_THREADS"] = "1"

from firedrake import *

# define mesh


# define function space, trial (u), and test (v) functions


# define forcing function f


# define bilinear and linear froms


# setup solution and solve, bcs not given since weakly enforced in our variational form


# write solution
