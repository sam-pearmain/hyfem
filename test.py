import abc

from typing import List, Mapping

from hyfem.utils import todo

class Linear(abc.ABC):
    def a(self):
        """returns the bilinear form a(u, v)"""

    def L(self):
        """returns the linear form L(v)"""

class Nonlinear(abc.ABC):
    def F(self):
        """returns the residual F(u) = 0"""

class TimeDependent(abc.ABC):
    """we have to be linear in time"""

class Equation(abc.ABC):
    """This is our base, we assume closed"""
    @abc.abstractmethod
    def state_variable(self) -> str:
        """return 'u' for example"""

class UnclosedEquation(abc.ABC):
    @abc.abstractmethod
    def requires(self): ...

    @abc.abstractmethod
    def solves_for(self): ...

class System:
    _equations: List[Equation]
    # the forms need to know where to draw their functions from, this maps the variable name to the function space
    _spaces: Mapping[str, FunctionSpace]

    def __init__(self, *eqns: Equation):
        self._equations = list(eqns)
        self._validate_closure()

    def state_variables(self): ...

    def _validate_closure():
        todo()