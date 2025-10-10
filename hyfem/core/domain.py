import ufl
import numbers

from firedrake import FunctionSpace, SpatialCoordinate, FacetNormal
from firedrake.mesh import MeshGeometry, MeshTopology


class Domain(object):
    """The problem domain"""
    def __init__(
            self, 
            mesh: MeshGeometry, 
            family: str, 
            degree: int, 
            name: str | None = None
        ) -> None:
        if not _has_standard_topology_backend(mesh):
            raise TypeError(f"unsupported mesh topology of type: {type(mesh.topology)}")
        
        self.mesh = mesh
        self.family = family
        self.degree = degree
        
        if name is not None:
            self.name = name 
        elif not _has_firedrake_default_name(self.mesh):
            self.name = self.mesh.name
        else:
            self.name = "_default_domain_name_"

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