import abc
import ufl
import ufl.equation

from typing import Self

from hyfem.core.discretisation import Discretisation
from hyfem.core.eqns.traits import Linear, Nonlinear, Solvable
from hyfem.utils import *


class Equation(Solvable):
    _discretisation_method: Discretisation | None

    def __init__(self) -> None:
        super().__init__()
        self._discretisation_method = None

    @ClassMethod
    def using_discretisation_method(cls, method: str) -> Self:
        eqn = cls()
        eqn.set_discretisation_method(method)
        return eqn

    @Property
    def discretisation_method(self) -> Discretisation:
        if self._discretisation_method is None:
            raise RuntimeError(f"discretisation scheme not set for {type(self).__name__}")
        return self._discretisation_method

    @Property
    def default_discretisation_method(self) -> Discretisation:
        return self._default_discretisation_method_impl()

    def set_discretisation_method(self, method: str) -> None:
        self._discretisation_method = Discretisation.from_str(method)

    def _is_system_of_equations_impl(self) -> bool:
        return False
    
    @abc.abstractmethod
    def _default_discretisation_method_impl(self) -> Discretisation:
        ...

class LinearEquation(Linear, Equation):
    @Property
    def a(self) -> ufl.Form:
        return self._bilinear_form_dispatch()
    
    @Property
    def L(self) -> ufl.Form:
        return self._linear_functional_dispatch()

    @Property
    def bilinear_form(self) -> ufl.Form:
        return self._bilinear_form_dispatch()
    
    @Property
    def linear_functional(self) -> ufl.Form:
        return self._linear_functional_dispatch()

    def _ufl_form_impl(self) -> ufl.Form:
        return self.a - self.L
    
    def _ufl_equation_impl(self) -> ufl.equation.Equation:
        return self.a == self.L

    def _bilinear_form_dispatch(self) -> ufl.Form:
        method = f"_{self.discretisation_scheme.abbrev()}_bilinear_form_impl"
        return getattr(self, method)()

    def _linear_functional_dispatch(self) -> ufl.Form:
        method = f"_{self.discretisation_scheme.abbrev()}_linear_functional_impl"
        return getattr(self, method)()

class NonlinearEquation(Nonlinear, Equation):
    @Property
    def F(self) -> ufl.Form:
        return self._nonlinear_form_dispatch()

    @Property
    def nonlinear_form(self) -> ufl.Form:
        return self._nonlinear_form_dispatch()

    def _ufl_form_impl(self) -> ufl.Form:
        return self._nonlinear_form_dispatch()
    
    def _ufl_equation_impl(self) -> ufl.equation.Equation:
        return self.F == 0

    def _nonlinear_form_dispatch(self) -> ufl.Form:
        method = f"_{self.discretisation_scheme.abbrev()}_nonlinear_form_impl"
        return getattr(self, method)()