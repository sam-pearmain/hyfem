import abc
import ufl
import ufl.equation

from typing import List

from hyfem.core.spaces import Spaces
from hyfem.utils import Property


class Solvable(abc.ABC):
    """The abstract base class for any equation or system of equations"""
    _spaces: Spaces | None

    def __init__(self) -> None:
        super().__init__()
        self._spaces = None

    @Property
    def ufl_form(self) -> ufl.Form: 
        return self._ufl_form_impl()
    
    @Property
    def ufl_equation(self) -> ufl.equation.Equation:
        return self._ufl_equation_impl()

    @Property
    def unknowns(self) -> List[str]:
        return self._unknowns_impl()

    def is_system_of_equations(self) -> bool:
        return self._is_system_of_equations_impl()

    def assign_function_spaces(self, spaces: Spaces) -> None:
        supported_eqn, _ = spaces.defined_on_str()
        if not supported_eqn == type(self).__name__:
            raise ValueError(f"attempted to assign <Spaces ({supported_eqn})> object to {type(self).__name__}")
        self._spaces = spaces

    def has_function_spaces(self) -> bool:
        return self._spaces is not None

    def has_defined_function_spaces(self) -> bool:
        if not self.has_function_spaces():
            raise RuntimeError(f"{type(self).__name__} has not been assigned a {Spaces.__name__} object")
        return self._spaces.all_spaces_assigned()

    def _ensure_defined_function_spaces(self) -> None:
        if not self.has_defined_function_spaces():
            raise RuntimeError(f"{type(self).__name__} has undefined function spaces") 

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

class LinearSolveable(Solvable):
    @Property
    def a(self) -> ufl.Form:
        return self._bilinear_form_impl()
    
    @Property
    def L(self) -> ufl.Form:
        return self._linear_functional_impl()

    @Property
    def bilinear_form(self) -> ufl.Form:
        return self._bilinear_form_impl()
    
    @Property
    def linear_functional(self) -> ufl.Form:
        return self._linear_functional_impl()

    @abc.abstractmethod
    def _bilinear_form_impl(self) -> ufl.Form:
        ...

    @abc.abstractmethod
    def _linear_functional_impl(self) -> ufl.Form:
        ...

    def _ufl_form_impl(self) -> ufl.Form:
        return self.a - self.L
    
    def _ufl_equation_impl(self) -> ufl.equation.Equation:
        return self.a == self.L

class NonlinearSolveable(Solvable):
    @Property
    def F(self) -> ufl.Form:
        return self._nonlinear_form_impl()

    @Property
    def nonlinear_form(self) -> ufl.Form:
        return self._nonlinear_form_impl()

    def _ufl_form_impl(self) -> ufl.Form:
        return self._nonlinear_form_impl()
    
    def _ufl_equation_impl(self) -> ufl.equation.Equation:
        return self.F == 0

    @abc.abstractmethod
    def _nonlinear_form_impl(self) -> ufl.Form:
        ...

class TimeDependent(abc.ABC):
    """A mixin for solvables with a time dimension"""
    ...

class PsuedotimeMixin(abc.ABC):
    ...