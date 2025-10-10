from firedrake import *


mesh = RectangleMesh(201, 101, 2.0, 1.0)
t = 0.0
t_end = 2.0
dt = 0.02

V = FunctionSpace(mesh, "DG", 1)

u_t = Function(V, name = "sol")
x, y = SpatialCoordinate(mesh)
u_init = 2 / sqrt(2 * 3.14159) * exp(-(1 / 2) * pow(5 * (x - 0.2), 2))
project(u_init, u_t) # set up initial condition, u_0

u = TrialFunction(V)
v = TestFunction(V)

c = as_vector((0.5, 0.0)) # constant velocity vector

n = FacetNormal(mesh)
u_flux = (
    0.5 * (dot(c, n('+')) + abs(dot(c, n('+')))) * u('+') + 
    0.5 * (dot(c, n('+')) - abs(dot(c, n('+')))) * u('-')
)

F_spatial = -dot(c, grad(v)) * u * dx + dot(c, n('+')) * u_flux * jump(v) * dS
F = (inner((u - u_t) / dt, v) * dx) + F_spatial

a = lhs(F)
L = rhs(F)

u_t_plus_one = Function(V)

problem = LinearVariationalProblem(a, L, u_t_plus_one)
solver = LinearVariationalSolver(
    problem, 
    solver_parameters = {
        "ksp_type": "preonly", 
        "pc_type":  "bjacobi",
        "sub_pc_type": "ilu"
        }, 
)

outfile = VTKFile("advection.pvd")
outfile.write(u_t)

while t < t_end:
    t += dt
    print(f"time: {t}/{t_end}")

    solver.solve()

    u_t.assign(u_t_plus_one)

    outfile.write(u_t)

print("finished")