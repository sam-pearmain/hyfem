
from firedrake.mesh import (MeshGeometry)
from hyfem.core.pde import (PDE, System)
from hyfem.core.spaces import Spaces


class Domain:
    _mesh: MeshGeometry
    _spaces: Spaces
    _equation: PDE | System

    def __init__(
            self, 
            mesh: MeshGeometry, 
            eqn: PDE | System, 
            spaces: Spaces,
            name: str
        ):
        self._mesh = mesh
        self._spaces = Spaces()
        self._equation = None
        
        if not self._has_standard_topology_backend():
            raise self._unsupported_mesh_topology()
        
        if name is not None:
            self.name = name 
        elif not self._has_firedrake_default_name():
            self.name = self.mesh.name
        else:
            self.name = f"_default_{self._equation.__name__}_domain_"