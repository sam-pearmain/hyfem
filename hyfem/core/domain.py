from typing import Optional, TypeVar, Generic

from firedrake import Function
from firedrake.mesh import MeshGeometry
from hyfem.firedrake.constructors import facet_normal
from hyfem.core.spaces import Spaces
from hyfem.core.mesh import Mesh
from hyfem.utils import *

if type_checking():
    from hyfem.core.eqns import Solvable


E = TypeVar('E', bound = 'Solvable')
class Domain(Generic[E]):
    _name: Optional[String]
    _mesh: MeshGeometry
    _spaces: Spaces
    _equation: E

    def __init__(
            self, 
            mesh: MeshGeometry, 
            eqn: E, 
            name: str | None = None
        ) -> None:
        self.name = name if name else f"_{eqn}_domain_"
        self._mesh = Mesh(mesh)
        self._spaces = Spaces(eqn, self._mesh)
        self._equation = eqn
        self._equation.assign_spaces(self._spaces)
        self._equation.assign_mesh(self._mesh)

        
    @Property
    def spaces(self) -> Spaces[E]:
        return self._spaces

    @Property
    def equation(self) -> E:
        return self._equation
    
    @Property
    def solution(self) -> Function:
        raise NotImplementedError("hands off")
    
def tests():
    pass

if __name__ == "__main__":
    tests()