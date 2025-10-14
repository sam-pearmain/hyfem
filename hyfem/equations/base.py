import abc

from typing import List

from hyfem.core.discretisation import DiscontinuousGalerkinMixin, DiscretisationScheme, ContinuousGalerkinMixin
from hyfem.utils import *


class Equation(abc.ABC):
    _discretisation: str

    def __init__(self, discretisation: str = "continuous galerkin"):
        if discretisation:
            self._discretisation = discretisation

    @Property
    def state_variables(self) -> List[str]: 
        """The state variables of the equation, for example [rho, rho_u, E]"""
        return self._state_variables_impl()

    @Property
    def auxiliary_variables(self) -> List[str]: 
        """The auxiliary variables of the equation, for example [k, ε]"""
        return self._auxiliary_variables_impl()

    @Property
    def variables(self) -> List[str]:
        """All the variables of the equation"""
        return self._state_variables_impl() + self._auxiliary_variables_impl()

    @Property
    def n_state_variables(self) -> int: 
        return len(self.state_variables())

    @Property
    def n_auxiliary_variables(self) -> int: 
        return len(self._auxiliary_variables_impl()) if self._auxiliary_variables_impl() else 0 

    @Property
    def n_variables(self) -> int:
        return len(self.variables)

    @ClassMethod
    def supported_discretisation_schemes(cls) -> List[DiscretisationScheme]:
        return [
            scheme for c in cls.__mro__
            if (scheme := DiscretisationScheme.from_cls(c)) is not None 
        ]


    @abc.abstractmethod
    def _state_variables_impl(self) -> List[str]: ...

    @abc.abstractmethod
    def _auxiliary_variables_impl(self) -> List[str] | None: ...


def tests():
    class Advection(Equation):
        """∂phi/∂t + div(phi u) = 0"""
        def _state_variables_impl(self):
            return ["phi"]
        
        def _auxiliary_variables_impl(self):
            return []
        
        def _cg_residual_impl(*args):
            return super()._cg_residual_impl()
        
        def _dg_residual_impl(*args):
            return super()._dg_residual_impl()
        
    class Poisson(Equation):
        def _state_variables_impl(self):
            return 
        
        def _auxiliary_variables_impl(self):
            return super()._auxiliary_variables_impl()
        
    eqn = Advection()
    print(eqn.state_variables)
    print(eqn.variables)
    print(eqn.n_variables)
    print(eqn.supported_discretisation_schemes())

    eqn2 = Poisson("mixed")
        
if __name__ == "__main__":
    tests()