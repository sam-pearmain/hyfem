import abc

from typing import List

from hyfem.core.eqns.eqn import LinearEquation
from hyfem.core.eqns.traits import LinearSolveable, NonlinearSolveable
from hyfem.core.spaces import Spaces

from hyfem.utils import *

class CoupledEquations(abc.ABC):
    ...

class LinearCoupledEquations(LinearSolveable, CoupledEquations):
    _spaces: Spaces | None
    _equations: List[LinearEquation]

    def __init__(self, eqns: List[LinearEquation]):
        super().__init__()
        self._equations = eqns

    def _unknowns_impl(self) -> List[str]:
        unknowns = []
        for eqn in self._equations:
            if eqn.unknown not in unknowns: unknowns.append(eqn.unknown)

    def is_system(self) -> bool: return True
    
class NonlinearCoupledEquations(NonlinearSolveable, CoupledEquations):
    ...