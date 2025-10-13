import ufl
import numbers

from typing import Any, Mapping, TypeVar, Generic, List
from firedrake import (
    SpatialCoordinate, FacetNormal, FunctionSpace, VectorFunctionSpace, BaseFunctionSpace, 
    MixedFunctionSpace
)
from firedrake.mesh import MeshGeometry, MeshTopology

from hyfem.equations.base import Equation
from hyfem.utils import *

E = TypeVar('E', bound = Equation)
class Domain(Generic[E], object):
    """The problem domain"""
    _mesh: MeshGeometry
    _spaces: Mapping[str, BaseFunctionSpace | None] # should probably be its own class
    _equation: Equation

    def __init__(
            self, 
            mesh: MeshGeometry, 
            equation: E, 
            name: str | None = None
        ) -> None:
        self._mesh = mesh
        self._spaces = {v: None for v in equation.variables()}
        self._equation = equation
        
        if not self._has_standard_topology_backend():
            raise TypeError(
                f"unsupported mesh topology of type: {type(self._mesh.topology)}"
            )
        
        if name is not None:
            self.name = name 
        elif not self._has_firedrake_default_name():
            self.name = self.mesh.name
        else:
            self.name = f"_default_{self._equation.__name__}_domain_"

    def assign_function_space(
            self, 
            family: str, 
            degree: int, 
            var_name: str, 
            vector_valued: bool = False, 
        ) -> None:
        """
        Assigns a variable from the domain's equation set to a given scalar-valued 
        function space
        """
        if (var_name not in self._equation.state_variables() and 
            var_name not in self._equation.auxiliary_variables()):
            raise self._var_not_found(var_name)
        
        if vector_valued:
            V = VectorFunctionSpace(self._mesh, family, degree, name = f"V_{var_name}")
        else:
            V = FunctionSpace(self._mesh, family, degree, name = f"V_{var_name}")

        if self._spaces.get(var_name) is not None:
            raise self._space_already_defined()
        
        self._spaces[var_name] = V
    
    def function_space(self, var_name: str | None = None) -> Any:
        """
        Returns the function space assigned to a given variable, if no variable is given
        then the entire mixed function space is returned
        """
        if var_name is not None:
            if var_name not in self._spaces:
                raise self._var_not_found(var_name)
            
            space = self._spaces[var_name]

            if space is None:
                raise self._space_not_defined(var_name)
            
            return space
        else:
            spaces = self._spaces.values()

            if any(spaces) == None:
                raise self._space_not_defined()

            return MixedFunctionSpace(spaces)

    @property
    def spatial_coordinates(self) -> SpatialCoordinate:
        return SpatialCoordinate(self._mesh)

    @property
    def facet_normal(self) -> ufl.FacetNormal:
        return FacetNormal(self._mesh)
    
    @property
    def topological_dimensions(self) -> numbers.Integral:
        return self._mesh.topological_dimension()
    
    @property
    def geometric_dimensions(self) -> numbers.Integral:
        return self._mesh.geometric_dimension()
    
    @check
    def _check_all_spaces_assigned(self) -> bool:
        if any(space is None for space in self._spaces.values()):
            return False
        return True
    
    @check
    def _has_standard_topology_backend(self) -> bool:
        return type(self._mesh.topology) is MeshTopology

    @check
    def _has_firedrake_default_name(self) -> bool:
        return self._mesh.name == "firedrake_default"

    @error
    def _var_not_found(self, var_name: str) -> ValueError:
        return ValueError(
            f"given variable: {var_name} not in {self._equation.__name__}'s\n" +
            f"state variables: {self._equation.state_variables()} \n or \n " +
            f"auxiliary variables: {self._equation.auxiliary_variables()}"
        )
    
    @error
    def _space_not_defined(self, *vars: str) -> ValueError:
        label = "variable" if len(vars) == 1 else "variables"
        names = ", ".join(vars)
        return ValueError(f"function space not defined for {label}: {names}")
    
    @error
    def _space_already_defined(self, *vars: str) -> ValueError:
        label = "variable" if len(vars) == 1 else "variables"
        details = ", ".join([f"{v} ∈ {self._spaces[v].name}" for v in vars])
        return ValueError(f"function space already defined for {label}: {details}")


def tests() -> None:
    from firedrake import UnitSquareMesh, ExtrudedMesh, VertexOnlyMesh

    mesh_standard = UnitSquareMesh(4, 4, name = "standard")
    mesh_extruded = ExtrudedMesh(UnitSquareMesh(4, 4), layers = 10, name = "extruded")
    mesh_vtx_only = VertexOnlyMesh(UnitSquareMesh(4, 4), [[0.5, 0.5]], name = "vtx_only")
    meshes: List[MeshGeometry] = [mesh_standard, mesh_extruded, mesh_vtx_only]

    filtered_meshes = [
        m for m in meshes 
        if (type(m.topology) is MeshTopology)
    ]

    assert [m.name for m in filtered_meshes] == ["standard"] 

if __name__ == "__main__":
    tests()