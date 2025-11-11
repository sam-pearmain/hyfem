import ufl
import numbers

from typing import Any

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

    def __getattr__(self, name) -> Any:
        return getattr(self._mesh, name)

    def __str__(self) -> String:
        return f"{type(self).__name__}".lower()

    @Property
    def ufl_domain(self) -> MeshGeometry:
        """Access the underlying ufl representation of the mesh"""
        return self._mesh

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
        return SpatialCoordinate(self.ufl_domain)

    @Property
    def facet_normal(self) -> ufl.FacetNormal:
        return facet_normal(self.ufl_domain)
    
    @Property
    def topological_dimensions(self) -> numbers.Integral:
        return self.ufl_domain.topological_dimension()
    
    @Property
    def geometric_dimensions(self) -> numbers.Integral:
        return self.ufl_domain.geometric_dimension()

    def _has_standard_topology_backend(self) -> bool:
        return type(self.ufl_domain.topology) is MeshTopology

    def _has_firedrake_default_name(self) -> bool:
        return self.ufl_domain.name == "firedrake_default"