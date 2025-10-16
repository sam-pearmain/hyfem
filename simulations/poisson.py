import os
os.environ["OMP_NUM_THREADS"] = "1"

from firedrake import *

n = 128
mesh = UnitSquareMesh(n, n)

V = FunctionSpace(mesh, "CG", 1)

u = TrialFunction(V)
v = TestFunction(V)

x, y = SpatialCoordinate(mesh)
f = Function(V).interpolate(sin(pi * x) * sin(pi * y))

a = inner(grad(v), grad(u)) * dx
L = f * v * dx

# bc = DirichletBC(V, Constant(0.0), "on_boundary")

u_sol = Function(V, name = "solution")
solve(a == L, u_sol)

outputfile = VTKFile("poisson.pvd")
outputfile.write(u_sol)
