import abc

from typing import List, Self, Type

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
        
    class Poisson(Equation):
        def _state_variables_impl(self):
            return 
        
        def _auxiliary_variables_impl(self):
            return super()._auxiliary_variables_impl()
        
    eqn = Advection()
    print(eqn.state_variables)
    print(eqn.variables)
    print(eqn.n_variables)
        
if __name__ == "__main__":
    tests()