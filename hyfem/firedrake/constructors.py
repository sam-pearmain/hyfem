from firedrake import (UnitIntervalMesh, UnitSquareMesh, UnitCubeMesh, 
                       FunctionSpace, VectorFunctionSpace, MixedFunctionSpace, 
                       TestFunction, TrialFunction, FacetNormal)

__all__ = [
    "unit_interval_mesh",
    "unit_square_mesh",
    "unit_cube_mesh",
    "function_space",
    "vector_function_space",
    "mixed_function_space",
    "test_function",
    "trial_function",
    "facet_normal", 
]

unit_interval_mesh = UnitIntervalMesh
unit_square_mesh = UnitSquareMesh
unit_cube_mesh = UnitCubeMesh
function_space = FunctionSpace
vector_function_space = VectorFunctionSpace
mixed_function_space = MixedFunctionSpace
test_function = TestFunction
trial_function = TrialFunction
facet_normal = FacetNormal
