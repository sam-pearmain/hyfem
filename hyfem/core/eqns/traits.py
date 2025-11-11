import abc
import ufl

from typing import List, Optional

from hyfem.utils import *

if type_checking():
    from hyfem.core.spaces import Spaces
    from hyfem.core.mesh import Mesh


class Solvable(abc.ABC):
    """The abstract base class for any equation or system of equations"""
    _spaces: Optional['Spaces']
    _mesh: Optional['Mesh']

    def __init__(self) -> None:
        super().__init__()
        self._spaces = None
        self._mesh = None

    def __str__(self) -> str:
        return f"{type(self).__name__}".lower()

    @Property
    def spaces(self) -> 'Spaces':
        if not self.has_spaces():
            raise AttributeError(f"spaces not initialised for {self}")
        return self._spaces
    
    @Property
    def mesh(self) -> 'Mesh':
        if not self.has_mesh():
            raise AttributeError(f"mesh not initialised for {self}")

    @Property
    def ufl_form(self) -> ufl.Form: 
        return self._ufl_form_impl()
    
    @Property
    def ufl_equation(self) -> ufl.equation.Equation:
        return self._ufl_equation_impl()

    @Property
    def unknowns(self) -> List[str]:
        return self._unknowns_impl()

    def is_solvable(self) -> bool:
        return True

    def is_system_of_equations(self) -> bool:
        return self._is_system_of_equations_impl()

    def assign_spaces(self, spaces: 'Spaces') -> None:
        defined_on_eqn, _ = spaces.defined_on()
        if not type(defined_on_eqn).__name__ == type(self).__name__:
            raise ValueError(f"attempted to assign <Spaces ({defined_on_eqn})> object to {self}")
        
        if self.has_spaces():
            raise AttributeError(f"spaces already set for {self}")

        self._spaces = spaces

    def assign_mesh(self, mesh: 'Mesh') -> None:
        if self.has_mesh():
            raise AttributeError(f"mesh already assigned for {self}")

        self._mesh = mesh

    def has_spaces(self) -> bool:
        return self._spaces is not None

    def has_mesh(self) -> bool:
        return self._mesh is not None

    @abc.abstractmethod
    def _ufl_form_impl(self) -> ufl.Form: 
        ...

    @abc.abstractmethod
    def _ufl_equation_impl(self) -> ufl.equation.Equation:
        ...

    @abc.abstractmethod
    def _unknowns_impl(self) -> List[str]:
        ...

    @abc.abstractmethod
    def _is_system_of_equations_impl(self) -> bool:
        ...

class LinearMixin(abc.ABC):
    """A mixin for linear equations"""
    def is_linear(self)    -> bool: return True
    def is_nonlinear(self) -> bool: return False

    @abc.abstractmethod
    def _bilinear_form_impl(self) -> ufl.Form:
        ...

    @abc.abstractmethod
    def _linear_functional_impl(self) -> ufl.Form:
        ...

class NonlinearMixin(abc.ABC):
    """A mixin for nonlinear equations"""
    def is_linear(self)    -> bool: return False
    def is_nonlinear(self) -> bool: return True

    @abc.abstractmethod
    def _nonlinear_form_impl(self) -> ufl.Form:
        ...

Linear = LinearMixin
Nonlinear = NonlinearMixin

class TimeMixin(abc.ABC):
    """A mixin for equations with a time dimension"""
    ...

class PsuedotimeMixin(abc.ABC):
    """A mixin for equatoins with a psuedo-time dimension"""
    ...

class SourceMixin(abc.ABC):
    """A mixin for equations with source terms"""
    @Property
    def f(self) -> ufl.Form:
        return self._f

    @abc.abstractmethod
    def _f_impl(self) -> ufl.Form:
        ...