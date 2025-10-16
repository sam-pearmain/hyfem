import abc
import ufl

from typing import List
from firedrake import (
    SpatialCoordinate, FacetNormal, FunctionSpace, VectorFunctionSpace, BaseFunctionSpace, 
    MixedFunctionSpace
)
from firedrake.mesh import MeshGeometry

from hyfem.core.discretisation import DiscretisationScheme
from hyfem.core.domain import Domain
from hyfem.core.spaces import Spaces
from hyfem.utils import *

class Equation(abc.ABC):
    _mesh: MeshGeometry
    _spaces: Spaces
    _discretisation: DiscretisationScheme

    def __init__(self, discretisation: str, mesh: MeshGeometry):
        self._discretisation = DiscretisationScheme.from_str(discretisation)

        if self._discretisation not in self.supported_discretisation_schemes():
            raise self._unsupported_discretisation_scheme()

    @Property
    def form(self, spaces):
        """Returns the ufl form of the equation"""
        return getattr(self, f"_{self._discretisation.abbrev()}_form_impl")(spaces)

    @ClassMethod
    def state_variables(cls) -> List[str]: 
        """The state variables of the equation, for example [rho, rho_u, E]"""
        return getattr(cls, f"_{cls._discretisation.abbrev()}_state_variables_impl")()

    @Property
    def n_state_variables(cls) -> int: 
        return len(cls.state_variables)

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
