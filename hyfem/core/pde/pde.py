import abc
import ufl

from hyfem.core.pde.system import Solvable
from hyfem.core.spaces import Spaces
from hyfem.utils import *

class Equation(Solvable):
    """
    A closed-form equation with a single unknown variable
    """
    _spaces: Spaces | None

    def __init__(self) -> None:
        super().__init__()

    @Property
    def unknowns(self) -> str: 
        return self._unknown_impl()

    @Property
    def F(self) -> ufl.Form:
        """The spatial residual, F(u)"""
        return self._spatial_residual_impl()
    
    @Property
    def spatial_residual(self) -> ufl.Form:
        """The spatial residual, F(u)"""
        return self._spatial_residual_impl()
        

    
    def is_system(self) -> bool: return False

    abc.abstractmethod
    def _unknown_impl(self) -> str: ...

    abc.abstractmethod
    def _spatial_residual_impl(self) -> ufl.Form: ...