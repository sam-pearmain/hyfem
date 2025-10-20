import abc
import warnings

from typing import List, Type, Self

from hyfem.core.pde import PDE
from hyfem.core.pde.traits import Solvable, Unclosed
from hyfem.core.spaces import Spaces

from hyfem.utils import *

class System(Solvable):
    _spaces: Spaces
    _equations: List[PDE] | None = None

    def __init__(self, eqns: List[PDE]):
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

    @Property
    def unknowns(self) -> List[str]:
        unknowns = []
        for eqn in self._equations:
            if eqn.unknown not in unknowns: unknowns.append(eqn.unknown)

    def is_system(self) -> bool: return True
    