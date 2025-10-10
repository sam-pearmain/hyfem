import ufl
import numbers

from typing import Self, Tuple, Type, List
from firedrake import SpatialCoordinate, FacetNormal
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
            self.name = "_default_domain_name_"

    def assign_variable(self, var_name: str, family: str, degree: int) -> None:
        """Assigns a variable from the domain's equation set to a given function space"""
        if (var_name not in self._equation.state_variables() and 
            var_name not in self._equation.auxiliary_variables()):
            raise ValueError(
                f"given variable {var_name} not in {self._equation.__name__}'s\n" +
                f"state variables: {self._equation.state_variables()} \nor \n " +
                f"auxiliary variables: {self._equation.auxiliary_variables()}"
            )
        
        
            

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