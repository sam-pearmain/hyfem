import ufl
import numbers

from firedrake import SpatialCoordinate
from firedrake.mesh import MeshGeometry, MeshTopology
from hyfem.firedrake import facet_normal
from hyfem.utils import *

class Mesh:
    _mesh: MeshGeometry
    _name: str | None = None

    def __init__(
            self, 
            mesh: MeshGeometry, 
        ):
        self._mesh = mesh
        self._name = mesh.name

        if not self._has_standard_topology_backend():
            raise TypeError(f"unsupported mesh topology: {type(self._mesh.topology)}")

    @Property
    def X(self) -> SpatialCoordinate:
        return self.spatial_coordinates

    @Property
    def x(self) -> SpatialCoordinate:
        x, _, _ = self.X
        return x

    @Property
    def y(self) -> SpatialCoordinate:
        if self.topological_dimensions < 2:
            raise RuntimeError(
                f"tried to get y component of a domain with only " +
                f"{self.topological_dimensions} dimensions"
            )
        _, y, _ = self.X
        return y

    @Property
    def z(self) -> SpatialCoordinate:
        if self.topological_dimensions < 3:
            raise RuntimeError(
                f"tried to get z component of a domain with only " +
                f"{self.topological_dimensions} dimensions"
            )
        _, _, z = self.X
        return z

    @Property
    def n(self) -> ufl.FacetNormal:
        return self.facet_normal

    @Property
    def spatial_coordinates(self) -> SpatialCoordinate:
        return SpatialCoordinate(self._mesh)

    @Property
    def facet_normal(self) -> ufl.FacetNormal:
        return facet_normal(self._mesh)
    
    @Property
    def topological_dimensions(self) -> numbers.Integral:
        return self._mesh.topological_dimension()
    
    @Property
    def geometric_dimensions(self) -> numbers.Integral:
        return self._mesh.geometric_dimension()

    def _has_standard_topology_backend(self) -> bool:
        return type(self._mesh.topology) is MeshTopology

    def _has_firedrake_default_name(self) -> bool:
        return self._mesh.name == "firedrake_default"