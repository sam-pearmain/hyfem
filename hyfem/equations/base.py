import abc
import ufl

from typing import List

from hyfem.core.discretisation import DiscretisationScheme
from hyfem.core.domain import Domain
from hyfem.utils import *

# some equation descriptors for later use
class Linear:
    def is_linear(self) -> bool:    return True
    def is_nonlinear(self) -> bool: return False

class Nonlinear:
    def is_linear(self) -> bool:    return False
    def is_nonlinear(self) -> bool: return True

class SourceMixin(abc.ABC):
    _supported_sources: List[ufl.Form]

    def __init__(self, *args):
        super().__init__(*args)
        self._supported_sources = self._supported_sources_impl()

    @Property
    def f(self, *args) -> ufl.Form:
        self._f_impl(*args)
    
    @abc.abstractmethod
    def _f_impl(self, *args) -> ufl.Form: ...

    @abc.abstractmethod
    def _supported_sources_impl(self) -> List[Callable[..., ufl.Form]]: ...


class LinearSourceMixin(SourceMixin, Linear):
    @Property
    def f(self, x) -> ufl.Form:
        self._f_impl(x)
    
    @abc.abstractmethod()
    def _f_impl(self, x) -> ufl.Form: ...

class NonlinearSourceMixin(SourceMixin, Nonlinear):
    @Property
    def f(self, u, x) -> ufl.Form:
        self._f_impl(u, x)
    
    @abc.abstractmethod
    def _f_impl(self, u, x) -> ufl.Form: ...

class TimeDependent: ...

class Equation(abc.ABC):
    _discretisation: DiscretisationScheme

    def __init__(self, discretisation: str, domain: Domain):
        self._discretisation = DiscretisationScheme.from_str(discretisation)

        if self._discretisation not in self.supported_discretisation_schemes():
            raise self._unsupported_discretisation_scheme()

    @Property
    def form(self, spaces):
        """Returns the ufl form of the equation"""
        return getattr(self, f"_{self._discretisation.abbrev()}_form_impl")(spaces)

    @Property
    def state_variables(self) -> List[str]: 
        """The state variables of the equation, for example [rho, rho_u, E]"""
        return getattr(self, f"_{self._discretisation.abbrev()}_state_variables_impl")()

    @Property
    def auxiliary_variables(self) -> List[str] | None: 
        """The auxiliary variables of the equation, for example [k, ε]"""
        return getattr(self, f"_{self._discretisation.abbrev()}_auxiliary_variables_impl")()

    @Property
    def variables(self) -> List[str]:
        """All the variables of the equation"""
        return self.state_variables + self.auxiliary_variables

    @Property
    def n_state_variables(self) -> int: 
        return len(self.state_variables)

    @Property
    def n_auxiliary_variables(self) -> int: 
        return len(self.auxiliary_variables) if self.auxiliary_variables else 0 

    @Property
    def n_variables(self) -> int:
        return len(self.variables)

    @ClassMethod
    def supported_discretisation_schemes(cls) -> List[DiscretisationScheme]:
        return [
            scheme for c in cls.__mro__
            if (scheme := DiscretisationScheme.from_cls(c)) is not None 
        ]
    
    @error
    def _unsupported_discretisation_scheme(self, scheme: str) -> Exception:
        return ValueError(
            f"unsupported discretisation scheme: {scheme} " +
            f"for {type(self.__name__)}"
        )


def tests():
    class Advection(Equation):
        """∂phi/∂t + div(phi u) = 0"""
        def _state_variables_impl(self):
            return ["phi"]
        
        def _auxiliary_variables_impl(self):
            return []
        
        def _cg_residual_impl(*args):
            return super()._cg_residual_impl()
        
        def _dg_residual_impl(*args):
            return super()._dg_residual_impl()
        
    class Poisson(Equation):
        def _state_variables_impl(self):
            return 
        
        def _auxiliary_variables_impl(self):
            return super()._auxiliary_variables_impl()
        
    eqn = Advection()
    print(eqn.state_variables)
    print(eqn.variables)
    print(eqn.n_variables)
    print(eqn.supported_discretisation_schemes())

    eqn2 = Poisson("mixed")
        
if __name__ == "__main__":
    tests()