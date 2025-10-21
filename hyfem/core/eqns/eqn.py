import abc
import ufl

from hyfem.core.spaces import Spaces
from hyfem.core.eqns.traits import LinearSolveable, NonlinearSolveable
from hyfem.utils import *

class LinearEquation(LinearSolveable):
    """
    A closed-form equation with a single unknown variable
    """
    _spaces: Spaces | None

    def __init__(self) -> None:
        super().__init__()

    @Property
    def unknowns(self) -> str: 
        return self._unknown_impl()

    abc.abstractmethod
    def _unknown_impl(self) -> str: ...

    abc.abstractmethod
    def _spatial_residual_impl(self) -> ufl.Form: ...

class NonlinearEquation(NonlinearSolveable):
    ...