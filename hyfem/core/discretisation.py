from enum import Enum, auto
from typing import Self, Type

from hyfem.utils import *

class DiscretisationScheme(Enum):
    """The enum of all discretisation schemes"""
    ContinuousGalerkin    = auto()
    DiscontinuousGalerkin = auto()

    @ClassMethod
    def from_str(cls, s: str) -> Self:
        match s.lower():
            case "continuous galerkin" | "cg":         return cls.ContinuousGalerkin
            case "discontinuous galerkin" | "dg":      return cls.DiscontinuousGalerkin
            case _: raise ValueError(f"unknown discretisation scheme: {s}")

    def __str__(self) -> str:
        match self:
            case self.ContinuousGalerkin:    return "continous galerkin"
            case self.DiscontinuousGalerkin: return "discontinuous galerkin"

    def __repr__(self) -> str: return str(self)

    def abbrev(self) -> str:
        match self:
            case self.ContinuousGalerkin:    return "cg"
            case self.DiscontinuousGalerkin: return "dg"