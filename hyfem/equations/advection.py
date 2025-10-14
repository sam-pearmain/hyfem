from typing import Mapping
from firedrake import *

from hyfem.core.discretisation import ContinuousGalerkinMixin
from hyfem.equations.base import Equation


class Advection(Equation, ContinuousGalerkinMixin):
    def _cg_state_variables_impl(self):
        return ["phi"]
    
    def _cg_auxiliary_variables_impl(self):
        return None
    
    def _cg_form_impl(self, spaces: Mapping[str, BaseFunctionSpace]):
        V = spaces[self._cg_state_variables_impl()[0]]
        
        u = TrialFunction(V)
        v = TestFunction(V)