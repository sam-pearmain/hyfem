import abc
import warnings

from typing import List, Type, Self

from hyfem.core.eqns.eqn import LinearEquation
from hyfem.core.eqns.traits import LinearSolveable, NonlinearSolveable
from hyfem.core.spaces import Spaces

from hyfem.utils import *

class LinearCoupledEquations(LinearSolveable):
    _spaces: Spaces | None
    _equations: List[LinearEquation]

    def __init__(self, eqns: List[LinearEquation]):
        super().__init__()
        self._equations = eqns

    @ClassMethod
    def of_equations(cls, eqns: List[PDE]) -> Self:
        for eqn in eqns:
            if not isinstance(eqn, PDE):
                raise TypeError("must pass list of pde instances")
        
            if Unclosed not in type(eqn).mro():
                warnings.warn("closed pde passed into system of equations, please check")
        
        return cls(eqns)

    def _unknowns_impl(self) -> List[str]:
        unknowns = []
        for eqn in self._equations:
            if eqn.unknown not in unknowns: unknowns.append(eqn.unknown)

    def is_system(self) -> bool: return True
    