from typing import Mapping, TypeVar

from firedrake import MeshGeometry
from firedrake.functionspaceimpl import FunctionSpace
from firedrake.functionspace import (FunctionSpace as function_space, 
                                     VectorFunctionSpace as vector_function_space) 

from hyfem.core.domain import Equation

E = TypeVar('E', bound = Equation)
class Spaces:
    _spaces: Mapping[str, FunctionSpace]
    
    def __init__(self) -> None:
        self._spaces = {}

    def assign_function_space(
            self, 
            eqn: Equation,
            mesh: MeshGeometry,
            family: str, 
            degree: int, 
            var_name: str,
            vector_valued: bool = False, 
        ) -> None:
        """
        Assigns a variable from the domain's equation set to a given scalar-valued 
        function space
        """
        if var_name not in eqn.state_variables():
            raise self._var_not_found(var_name)
        
        if vector_valued:
            V = vector_function_space(self._mesh, family, degree, name = f"V_{var_name}")
        else:
            V = function_space(self._mesh, family, degree, name = f"V_{var_name}")

        if self._spaces.get(var_name) is not None:
            raise self._space_already_defined()