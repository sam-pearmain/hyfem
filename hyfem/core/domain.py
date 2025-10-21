
from typing import TypeVar, Generic

from firedrake.mesh import (MeshGeometry, MeshTopology)
from hyfem.core.eqns import Solvable
from hyfem.core.spaces import Spaces
from hyfem.utils import *

E = TypeVar('E', bound = Solvable)
class Domain(Generic[E]):
    _mesh: MeshGeometry
    _spaces: Spaces
    _eqn: E

    def __init__(
            self, 
            mesh: MeshGeometry, 
            eqn: E, 
            name: str | None = None
        ) -> None:
        self._mesh = mesh
        self._spaces = Spaces(eqn, mesh)
        self._eqn = eqn
        self._eqn.assign_function_spaces(self._spaces)
        
        if not self._has_standard_topology_backend():
            raise TypeError(f"unsupported mesh topology: {type(self._mesh.topology)}")
        
        if name is not None:
            self.name = name 
        elif not self._has_firedrake_default_name():
            self.name = self.mesh.name
        else:
            self.name = f"_default_{self._eqn.__name__}_domain_"

    def _has_standard_topology_backend(self) -> bool:
        return type(self._mesh.topology) is MeshTopology

    def _has_firedrake_default_name(self) -> bool:
        return self._mesh.name == "firedrake_default"
   
    
def tests():
    pass

if __name__ == "__main__":
    tests()