import ufl
import numbers

from typing import Self, Tuple, Type, List
from firedrake import SpatialCoordinate, FacetNormal, FunctionSpace, VectorFunctionSpace
from firedrake.mesh import MeshGeometry, MeshTopology
from hyfem.equations.base import Equation
from hyfem.core.spaces import Spaces

class Domain(object):
    """The problem domain"""

    def __init__(
            self, 
            mesh: MeshGeometry, 
            equation: Type[Equation], 
            name: str | None = None
        ) -> None:
        if not _has_standard_topology_backend(mesh):
            raise TypeError(f"unsupported mesh topology of type: {type(mesh.topology)}")
        
        self._mesh = mesh
        self._spaces = []
        self._equation = equation
        
        if name is not None:
            self.name = name 
        elif not _has_firedrake_default_name(self.mesh):
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
        """Assigns a variable from the domain's equation set to a given scalar-valued function space"""
        if (var_name not in self._equation.state_variables() and 
            var_name not in self._equation.auxiliary_variables()):
            raise ValueError(
                f"given variable: {var_name} not in {self._equation.__name__}'s\n" +
                f"state variables: {self._equation.state_variables()} \nor \n " +
                f"auxiliary variables: {self._equation.auxiliary_variables()}"
            )
        
        if vector_valued:
            V = VectorFunctionSpace(self._mesh, family, degree, name = f"V_{var_name}")
        else:
            V = FunctionSpace(self._mesh, family, degree, name = f"V_{var_name}")

        if V in self._spaces:
            raise ValueError(f"function space {V} already defined")

        self._spaces.append(V)
    
    @property
    def spatial_coordinates(self) -> SpatialCoordinate:
        return SpatialCoordinate(self.mesh)

    @property
    def facet_normal(self) -> ufl.FacetNormal:
        return FacetNormal(self.mesh)
    
    @property
    def topological_dimensions(self) -> numbers.Integral:
        return self.mesh.topological_dimension()
    
    @property
    def geometric_dimensions(self) -> numbers.Integral:
        return self.mesh.geometric_dimension()
    
def _has_standard_topology_backend(mesh: MeshGeometry) -> bool:
    return type(mesh.topology) is MeshTopology

def _has_firedrake_default_name(mesh: MeshGeometry) -> bool:
    return mesh.name == "firedrake_default"


def tests() -> None:
    from firedrake import UnitSquareMesh, ExtrudedMesh, VertexOnlyMesh

    mesh_standard = UnitSquareMesh(4, 4, name = "standard")
    mesh_extruded = ExtrudedMesh(UnitSquareMesh(4, 4), layers = 10, name = "extruded")
    mesh_vtx_only = VertexOnlyMesh(UnitSquareMesh(4, 4), [[0.5, 0.5]], name = "vtx_only")
    meshes = [mesh_standard, mesh_extruded, mesh_vtx_only]

    filtered_meshes = [
        m for m in meshes 
        if _has_standard_topology_backend(m)
    ]

    assert [m.name for m in filtered_meshes] == ["standard"] 

if __name__ == "__main__":
    tests()