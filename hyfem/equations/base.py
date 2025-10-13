import abc
from typing import List
from firedrake import *

class Equation(abc.ABC):
    @classmethod
    @abc.abstractmethod
    def state_variables(cls) -> List[str]: 
        """The state variables of the equation, for example [rho, rho_u, E]"""
        ...

    @classmethod
    @abc.abstractmethod
    def auxiliary_variables(cls) -> List[str] | None: 
        """The auxiliary variables of the equation, for example [k, ε]"""
        ...

    @classmethod
    def variables(cls) -> List[str]:
        """All the variables of the equation"""
        return cls.state_variables() + cls.auxiliary_variables()

    @classmethod
    def n_state_variables(cls) -> int: return len(cls.state_variables())

    @classmethod
    def n_auxiliary_variables(cls) -> int: 
        if not cls.auxiliary_variables():
            return 0
        return len(cls.auxiliary_variables())


def tests():
    class Advection(Equation):
        """∂phi/∂t + div(phi u) = 0"""
        _form = ""

        def state_variables(cls) -> List[str]:
            return ["u"]
        
        def auxiliary_variables(cls) -> List[str] | None:
            return None
        
        def residual(cls):

        
        

if __name__ == "__main__":
    tests()