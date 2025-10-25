import abc
from re import T
import ufl
import ufl.equation

from typing import List

from hyfem.core.domain import Domain
from hyfem.utils import *


class Solvable(abc.ABC):
    """The abstract base class for any equation or system of equations"""
    _domain: Domain | None

    def __init__(self) -> None:
        super().__init__()
        self._domain = None

    @Property
    def domain(self) -> Domain:
        if self._domain is None:
            raise AttributeError(f"tried to access domain of {type(self).__name__} but it is unset")
        return self._domain
    
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

    def assign_domain(self, domain: Domain) -> None:
        domain_equations = type(domain.equation).__name__
        if not domain_equations == type(self).__name__:
            raise ValueError(f"attempted to assign <Spaces ({domain_equations})> object to {type(self).__name__}")
        self._domain = domain

    def has_domain(self) -> bool:
        return self._domain is not None

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