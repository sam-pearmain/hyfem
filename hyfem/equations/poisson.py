from firedrake import *

from hyfem.core.discretisation import ContinuousGalerkinMixin
from hyfem.core.eqns.eqn import LinearEquation
from hyfem.core.eqns.traits import SourceMixin
from hyfem.firedrake import unit_square_mesh
from hyfem.utils import *


class Poisson(LinearEquation, ContinuousGalerkinMixin, SourceMixin):
    def __init__(self):
        super().__init__()
    
    def _unknowns_impl(self):
        return ['U']
    
    def _bilinear_form_impl(self):
        u = self.spaces.get_trial_function('U')
        v = self.spaces.get_test_function('U')
        return inner(grad(u), grad(v)) * dx

    def _linear_functional_impl(self):
        v = self.spaces.get_test_function('U')
        return inner(self.f, v) * dx
    
    def _f_impl(self):
        V = self.spaces.get_function_space('U')
        x, y = self.mesh.spatial_coordinates
        return Function(V).interpolate(x + y)


def tests():
    from hyfem.core.domain import Domain
    from firedrake import DirichletBC, Constant

    mesh = unit_square_mesh(10, 10)
    equation = Poisson()
    domain = Domain(mesh, equation, "poisson")
    domain.spaces.assign_function_space('U', "CG", 1)
    solution = Function(domain.spaces.get_function_space('U'))
    bcs = DirichletBC(domain.spaces.get_function_space('U'), Constant(0), "on_boundary")

    problem = LinearVariationalProblem(
        domain.equation.a, 
        domain.equation.L, 
        solution, 
        bcs = bcs
    )
    solver = LinearVariationalSolver(problem)
    solver.solve()

    outfile = VTKFile("out.pvd")
    outfile.write(solution)

def test2():
    mesh = unit_square_mesh(10, 10)
    V = FunctionSpace(mesh, "CG", 2)
    u = TrialFunction(V)
    v = TestFunction(V)

    a = inner

if __name__ == "__main__":
    tests()