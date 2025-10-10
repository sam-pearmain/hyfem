import ufl
import numbers

from firedrake import FunctionSpace, SpatialCoordinate, FacetNormal
from firedrake.mesh import MeshGeometry, ExtrudedMeshTopology, VertexOnlyMeshTopology



class Domain(object):
    """
    Contains information about the problem domain. 
    """
    def __init__(
            self, 
            mesh: MeshGeometry, 
            family: str, 
            degree: int, 
            name: str | None = None
        ) -> None:
        if _is_standard_topology_backend(mesh):
            self.mesh = mesh
        else: 
            raise TypeError(f"unsupported mesh topology of type: {type(mesh.topology)}")
        
        self.family = family
        self.degree = degree
        
        if name is not None:
            self.name = name 
        elif self.mesh.name is not "firedrake_default":
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

    

def _is_standard_topology_backend(mesh: MeshGeometry) -> bool:
    """Checks whether a Firedrake mesh uses the standard MeshTopology backend"""
    topology = mesh.topology
    return not isinstance(topology, (ExtrudedMeshTopology, VertexOnlyMeshTopology))

def tests() -> None:
    from firedrake import UnitSquareMesh, ExtrudedMesh, VertexOnlyMesh

    mesh_standard = UnitSquareMesh(4, 4)
    mesh_extruded = ExtrudedMesh(UnitSquareMesh(4, 4), layers = 10, name = "extruded")
    mesh_vtx_only = VertexOnlyMesh(UnitSquareMesh(4, 4), [[0.5, 0.5]], name = "vtx_only")
    meshes = [mesh_standard, mesh_extruded, mesh_vtx_only]

    filtered_meshes = [
        m for m in meshes 
        if _is_standard_topology_backend(m)
    ]

    print(mesh_standard.name)
    # assert [m.name for m in filtered_meshes] == ["standard"] 

if __name__ == "__main__":
    tests()