import abc

from typing import Generic, List

from hyfem.core.eqns import Solvable
from hyfem.core.eqns.eqn import Equation, LinearEquation, Nonlinear
from hyfem.core.eqns.traits import Linear, Nonlinear, Solvable 
from hyfem.utils import *

E = TypeVar('E', bound = Equation)
class CoupledEquations(Generic[E], Solvable):
    _equations: List[E]

    def _is_system_of_equations_impl(self) -> bool:
        return True

LE = TypeVar('LE', bound = LinearEquation)
class LinearCoupledEquations(Linear, CoupledEquations[LE]):
    _equations: List[LE]

    def __init__(self, eqns: List[LE]):
        super().__init__()
        self._equations = eqns

E = TypeVar('E', bound = Equation)
class NonlinearCoupledEquations(Nonlinear, CoupledEquations[E]):
    ...