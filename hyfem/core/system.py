import abc
import warnings

from typing import List, Type, Self

from hyfem.core.pde import (PDE, Unclosed)
from hyfem.utils import *

class System(abc.ABC):
    _equations: List[Type[PDE]] # we don't need instances

    @ClassMethod
    def of_equations(cls, eqns: List[PDE]) -> Self:
        for eqn in eqns:
            if not isinstance(eqn, type):
                raise TypeError("must pass equation types not instances")
        
            if Unclosed not in eqn.mro():
                warnings.warn("closed pde passed into system of equations, please check")
        
        cls._equations = eqns

    @ClassMethod
    def unknowns(cls) -> List[str]:
        unknowns = []
        for eqn in cls._equations:
            if eqn.unknown not in unknowns: unknowns.append(eqn.unknown)

    @ClassMethod
    def 
    