import os
os.environ["OMP_NUM_THREADS"] = "1"

from firedrake import *
from firedrake.adjoint import *

pyadjoint.Tape.progress_bar = ProgressBar

continue_annotation()

n = 30
mesh = UnitIntervalMesh(n)
timestep = Constant(1.0 / n)
steps = 10

x, = SpatialCoordinate(mesh)
V = FunctionSpace(mesh, "CG", 2)
ic = project(sin(2.0 * pi * x), V, name = "ic")

u_new = Function(V, name = "u_new")
u_old = Function(V, name = "u_old")
u_old.assign(ic)

v = TestFunction(V)
nu = Constant(1e-4)
F = (
    (u_new - u_old) / timestep * v
    + u_new * u_new.dx(0) * v
    + nu * u_new.dx(0) * v.dx(0)
) * dx
bc = DirichletBC(V, 0.0, "on_boundary")
problem = NonlinearVariationalProblem(F, u_new, bcs = bc)
solver = NonlinearVariationalSolver(problem)

J = assemble(ic * ic * dx)

for _ in range(steps):
    solver.solve()
    u_old.assign(u_new)
    J += assemble(u_new * u_new * dx)

J_hat = ReducedFunctional(J, Control(ic))

ic_new = project(sin(pi*x), V)
J_new = J_hat(ic_new)

dJ = J_hat.derivative()

pause_annotation()
print(J, dJ)
continue_annotation()