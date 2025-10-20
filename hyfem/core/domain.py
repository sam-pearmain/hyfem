
from firedrake.mesh import (MeshGeometry, MeshTopology)
from hyfem.core.pde import (PDE, System)
from hyfem.core.spaces import Spaces
from hyfem.utils import *


class Domain:
    _mesh: MeshGeometry
    _spaces: Spaces
    _equation: PDE | System

    def __init__(
            self, 
            mesh: MeshGeometry, 
            eqn: PDE | System, 
            name: str | None = None
        ) -> None:
        self._mesh = mesh
        self._spaces = Spaces(eqn, mesh)
        self._equation = eqn
        
        if not self._has_standard_topology_backend():
            raise TypeError(f"unsupported mesh topology: {type(self._mesh.topology)}")
        
        if name is not None:
            self.name = name 
        elif not self._has_firedrake_default_name():
            self.name = self.mesh.name
        else:
            self.name = f"_default_{self._equation.__name__}_domain_"

    def _has_standard_topology_backend(self) -> bool:
        return type(self._mesh.topology) is MeshTopology

    def _has_firedrake_default_name(self) -> bool:
        return self._mesh.name == "firedrake_default"