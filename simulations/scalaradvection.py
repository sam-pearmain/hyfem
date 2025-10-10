from typing import Any, Self
from firedrake import *


class ScalarAdvection:
    """
    ∂u / ∂t + ∇·(uc) = f
    where,
        u: scalar solution field R^d -> R
        c: velocity vector field R^d -> R^d
        f: forcing term, often 0
    """
    
    def __enter__(self) -> Self:
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        pass

    def set_domain(self, domain) -> Self:
        self.domain = domain
        self.spatial_dimensions = domain.geometric_dimension()
        return self

    def set_function_space(self, family, degree) -> Self:
        self._function_space = 
        self.V = FunctionSpace(self.domain, family, degree)


        self.u, self.v = TrialFunction(V), TestFunction(V)

    def set_velocity_field(self, expr: Any) -> Self:
        pass

    def build_ufl(self):
        if self.function


def tests():
    from scalaradvection import ScalarAdvection

    mesh = UnitSquareMesh(10, 10)

    with ScalarAdvection() as eqn:
        eqn.set_domain(mesh)
        eqn.set_function_space("CG", 2)
        eqn.set_velocity_field()

    print(type(eqn.spatial_dimensions))
        

if __name__ == "__main__":
    tests()