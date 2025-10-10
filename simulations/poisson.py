from firedrake import * # type: ignore

n = 128
mesh = UnitSquareMesh(n, n)

V = FunctionSpace(mesh, "CG", 1)

u = TrialFunction(V)
v = TestFunction(V)

x, y = SpatialCoordinate(mesh)
f = Function(V).interpolate(sin(pi * x) * sin(pi * y))

a = dot(grad(v), grad(u)) * dx
L = v * f * dx

bc = DirichletBC(V, Constant(0.0), "on_boundary")

u_sol = Function(V)
solve(a == L, u_sol, bcs = bc)

outputfile = VTKFile("poisson.pvd")
outputfile.write(u_sol)
