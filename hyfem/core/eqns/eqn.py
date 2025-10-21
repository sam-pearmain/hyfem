import abc
import ufl
import ufl.equation

from hyfem.core.discretisation import DiscretisationScheme
from hyfem.core.eqns.traits import Linear, Nonlinear, Solvable
from hyfem.utils import *


class Equation(Solvable):
    _discretisation_scheme: DiscretisationScheme | None

    def __init__(self) -> None:
        super().__init__()
        self._discretisation_scheme = None

    @Property
    def discretisation_scheme(self) -> DiscretisationScheme:
        if self._discretisation_scheme is None:
            raise RuntimeError(f"discretisation scheme not set for {type(self).__name__}")
        return self._discretisation_scheme

    def _is_system_of_equations_impl(self) -> bool:
        return False

class LinearEquation(Linear, Equation):
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

    def _ufl_form_impl(self) -> ufl.Form:
        return self.a - self.L
    
    def _ufl_equation_impl(self) -> ufl.equation.Equation:
        return self.a == self.L

    def _bilinear_form_impl(self) -> ufl.Form:
        method = f"_{self.discretisation_scheme.abbrev()}_bilinear_form_impl"
        return getattr(self, method)()

    def _linear_functional_impl(self) -> ufl.Form:
        method = f"_{self.discretisation_scheme.abbrev()}_linear_functional_impl"
        return getattr(self, method)()

class NonlinearEquation(Nonlinear, Equation):
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

    def _nonlinear_form_impl(self) -> ufl.Form:
        method = f"_{self.discretisation_scheme.abbrev()}_nonlinear_form_impl"
        return getattr(self, method)()