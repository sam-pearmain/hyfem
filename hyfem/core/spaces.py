from typing import List, Type

from hyfem.core.pde import (PDE, System)
from hyfem.utils import *

class Spaces:
    _vars: List[str]

    def __init__(self, eqn: PDE | System) -> None:
        match type(eqn):
            case PDE():    self._vars = [eqn.unknown]
            case System(): self._vars = eqn.unknowns
            case _: raise TypeError(
                f"eqn must be either of type: {PDE.__name__} or {System.__name__}" + 
                f", not {type(eqn).__name__}"
            )