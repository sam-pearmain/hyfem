from firedrake import *

from hyfem.core.discretisation import ContinuousGalerkinMixin
from hyfem.core.eqns.eqn import LinearEquation
from hyfem.core.eqns.traits import SourceMixin
from hyfem.firedrake import trial_function, test_function
from hyfem.utils import *

if type_checking():
    from hyfem.core.domain import Domain


class Poisson(LinearEquation, ContinuousGalerkinMixin, SourceMixin):
    def __init__(self):
        super().__init__()
    
    def _unknowns_impl(self):
        return ['U']
    
    def _bilinear_form_impl(self):
        u = self.domain.spaces.get_trial_function('U')
        v = self.domain.spaces.get_test_function('U')
        return inner(grad(u), grad(v)) * dx

    def _linear_functional_impl(self):
        v = self.domain.spaces.get_test_function('U')
        return inner(self.f, v) * dx
    
    def _f_impl(self):
        V = self.domain.spaces.get_function_space('U')
        x, y = self.domain.spatial_coordinates
        return Function(V).interpolate(x + y)

def tests():
    from hyfem.core.domain import Domain

    mesh = UnitSquareMesh(10, 10)
    equation = Poisson()
    domain = Domain(mesh, equation, "poisson")
    domain.spaces.assign_function_space('U', "CG", 1)

    problem = LinearVariationalProblem(domain.equation.a, domain.equation.L)
    solver = LinearVariationalSolver(problem)
    solver.solve()

if __name__ == "__main__":
    tests()