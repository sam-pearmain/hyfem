import ufl

from typing import Mapping
from firedrake import *

from hyfem.core.eqns.eqn import LinearEquation
from hyfem.utils import *

if type_checking():
    from hyfem.core.domain import Domain


class Poisson(LinearEquation):
    def __init__(self, discretisation, f: ufl.Form | None = None):
        super().__init__(discretisation)
    
    def _cg_state_variables_impl(self):
        return ["u"]
    
    def _cg_auxiliary_variables_impl(self):
        return None
    
    def _cg_form_impl(self, domain: 'Domain'):
        V = domain.solution_space()

        u = TrialFunction(V)
        v = TestFunction(V)
        n = domain.facet_normal()

        a = inner(grad(u), grad(v)) * dx
        L = inner(self.f * v) * dx

