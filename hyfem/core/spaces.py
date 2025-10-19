from typing import List, Mapping, Type

from firedrake import BaseFunctionSpace

from hyfem.core.pde import (PDE, System)
from hyfem.utils import *

class Spaces:
    _eqn: PDE | System
    _vars: List[str]
    _spaces: Mapping[str, BaseFunctionSpace | None]

    def __init__(self, eqn: PDE | System) -> None:
        self._eqn = eqn
        match type(self._eqn):
            case PDE():    self._vars = [self._eqn.unknown]
            case System(): self._vars = self._eqn.unknowns
            case _: raise TypeError(
                f"eqn must be either of type: {PDE.__name__} or {System.__name__}" + 
                f", not {type(eqn).__name__}"
            )
        self._spaces = {var: None for var in self._vars}

    @Property
    def function_spaces(self) -> Mapping[str, BaseFunctionSpace]:
        if not self._spaces.v:
            raise AttributeError(
                f"no function spaces assigned"
            )
        return self._spaces
    
    @Property
    def function_space(self, var: str) -> BaseFunctionSpace:
        if var not in self._spaces.keys():
            raise ValueError(
                f"{var} not in the {type(self._eqn).__name__} unknowns"
            )
        
        return self._spaces[var]
    
    def _validate_variable(self, var) -> None:
        if var not in self._spaces.keys():
            raise ValueError(
                f"{var} not in the {type(self._eqn).__name__} unknowns"
            )

    def _space_assigned(self, var) -> bool:
        if var not in self._spaces.keys():
            raise ValueError(
                f"{var} not in the {type(self._eqn).__name__} unknowns"
            )
        
        return self._spaces[var] is None

    def _all_spaces_assigned(self) -> bool:
        return all(space is not None for space in self._spaces.values())