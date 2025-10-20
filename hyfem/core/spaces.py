import warnings

from typing import List, Mapping, Tuple
from numpy.typing import ArrayLike

from firedrake import BaseFunctionSpace, MeshGeometry, Function
from firedrake.functionspaceimpl import MixedFunctionSpace
from firedrake.ufl_expr import Argument, Coargument
from ufl.argument import Coargument

from hyfem.core.pde import (PDE, System)
from hyfem.core.pde.traits import Solvable
from hyfem.firedrake import *
from hyfem.utils import *

class Spaces(object):
    _eqn: PDE | System
    _mesh: MeshGeometry
    _vars: List[str]
    _spaces: Mapping[str, BaseFunctionSpace | None]

    def __init__(self, eqn: PDE | System, mesh: MeshGeometry) -> None:
        self._eqn = eqn
        self._mesh = mesh

        match type(self._eqn):
            case PDE():    self._vars = [self._eqn.unknown]
            case System(): self._vars = self._eqn.unknowns
            case _: raise TypeError(
                f"eqn must be either of type: {PDE.__name__} or {System.__name__}" + 
                f", not {type(eqn).__name__}"
            )

        self._spaces = {var: None for var in self._vars}

    def assign_function_space(
            self, 
            var: str, 
            family: str, 
            degree: int,
            name: str | None = None, 
            vector_valued: bool = False, 
        ) -> None:
        """
        Assigns a function space to the given variable
        """
        if self._space_assigned(var):
            raise RuntimeError(f"function space already assigned to {var} : {self._spaces[var]}")

        name = name if name else f"V_{var}"

        if vector_valued:
            space = vector_function_space(self._mesh, family, degree, name = name)
        else:
            space = function_space(self._mesh, family, degree, name = name)

        self._spaces[var] = space

    def update_function_space(
            self, 
            var: str, 
            family: str, 
            degree: int, 
            name: str | None = None, 
            vector_valued: bool = False, 
        ) -> None:
        """
        Updates the given variable's function space
        """
        if not self._space_assigned(var):
            warnings.warn(f"update_function_space called on {var} but space was previously unassigned")
            self.assign_function_space(var, family, degree, name, vector_valued)

        name = self._spaces[var].label()

        if vector_valued:
            space = vector_function_space(self._mesh, family, degree, name = name)
        else:
            space = function_space(self._mesh, family, degree, name = name)

        # this might cause some problems if we change from CG to DG, for example
        self._spaces[var] = space

    def defined_on(self) -> Tuple[Solvable, MeshGeometry]:
        return tuple(self._eqn, self._mesh)

    def defined_on_str(self) -> Tuple[str, str]:
        return tuple(type(self._eqn).__name__, type(self._mesh).__name__)

    def get_mixed_function_space(self, *vars: str, name: str | None = None) -> MixedFunctionSpace:
        """Returns the mixed space for the given vars, if no vars are given returns entire mixed space"""
        if not self._eqn.is_system():
            raise RuntimeError(
                f"{type(self._eqn).__name__} has only one unknown and " +
                f"therefore has no mixed function space"
            )
        
        spaces = self.get_function_spaces(*vars)
        return mixed_function_space(spaces, name) 
    
    def get_function_space(self, var: str) -> BaseFunctionSpace:
        """Returns the function space assigned to the given variable"""
        self._validate_variable(var)
        if not self._space_assigned(var):
            raise AttributeError(f"{var} has no assigned function space")
        return self._spaces[var]

    def get_function_spaces(self, *vars: str) -> Tuple[BaseFunctionSpace, ...]:
        """Gets the function spaces for the given vars, if no vars given returns all spaces"""
        if not vars:
            if not self.all_spaces_assigned():
                raise AttributeError(f"not all function spaces assigned")
            return tuple(self._spaces[var] for var in self._vars)

        return tuple(self.get_function_space(var) for var in vars)
    
    def get_trial_function(self, var: str) -> Coargument | Argument:
        """Creates a trial function in the given var's function space"""
        self._validate_variable(var)
        return trial_function(self._spaces[var])
    
    def get_trial_functions(self, *vars: str) -> Tuple[Coargument, ...] | Tuple[Argument, ...]:
        """Creates a list of trial functions in the given vars' function spaces"""
        spaces = self.get_function_spaces(*vars)
        return tuple(trial_function(space) for space in spaces)
    
    def get_test_function(self, var: str) -> Coargument | Argument:
        """Creates a test function in the given var's function space"""
        self._validate_variable(var)
        return test_function(self._spaces[var])
    
    def get_test_functions(self, *vars: str) -> Tuple[Coargument, ...] | Tuple[Argument, ...]:
        """Creates a list of test functions in the given vars' function spaces"""
        spaces = self.get_function_spaces(*vars)
        return tuple(test_function(space) for space in spaces)

    def create_function(self, var: str, val: ArrayLike | None = None, name: str | None = None) -> Function:
        """Creates and returns a Firedrake Function for the given variable"""
        space = self.get_function_space(var)
        return Function(space, val = val, name = name)
    
    def create_functions(self, *vars: str) -> Tuple[Function]:
        """Creates and returns a tuple of Firedrake Functions in the given vars' function spaces"""
        spaces = self.get_function_spaces(*vars)
        return tuple(Function(space) for space in spaces)

    def all_spaces_assigned(self) -> bool:
        """Checks whether all variables have an assigned function space"""
        return all(space is not None for space in self._spaces.values())

    def _validate_variable(self, var: str) -> None:
        """Validates whether the given variable exists within the spaces"""
        if var not in self._spaces.keys():
            raise ValueError(
                f"{var} not in the {type(self._eqn).__name__} unknowns"
            )

    def _space_assigned(self, var: str) -> bool:
        """Checks whether the given var already has an assigned function space"""
        return self._spaces[var] is not None