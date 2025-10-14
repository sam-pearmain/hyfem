import abc
import ufl

from enum import Enum, auto
from typing import List, Self, Type

from hyfem.utils import *

if type_checking():
    from hyfem.equations.base import Equation


class ContinuousGalerkinMixin(abc.ABC):
    @abc.abstractmethod
    def _cg_state_variables_impl(self) -> List[str]: ...

    @abc.abstractmethod
    def _cg_auxiliary_variables_impl(self) -> List[str] | None: ...

    @abc.abstractmethod
    def _cg_form_impl(self, *args): ...

class DiscontinuousGalerkinMixin(abc.ABC):
    @abc.abstractmethod
    def _dg_state_variables_impl(self) -> List[str]: ...

    @abc.abstractmethod
    def _dg_auxiliary_variables_impl(self) -> List[str] | None: ...

    @abc.abstractmethod
    def _dg_form_impl(self, *args): ...

class DiscretisationScheme(Enum):
    """The enum of all discretisation schemes"""
    ContinuousGalerkin    = auto()
    DiscontinuousGalerkin = auto()
    # MixedFormulation                = auto()
    # HybridisedDiscontinuousGalerkin = auto()

    @ClassMethod
    def from_str(cls, s: str) -> Self:
        match s.lower():
            case "continuous galerkin" | "cg":         return cls.ContinuousGalerkin
            case "discontinuous galerkin" | "dg":      return cls.DiscontinuousGalerkin
            case _: raise ValueError(f"unknown discretisation scheme: {s}")

    @ClassMethod
    def from_cls(cls, c: Type) -> Self | None:
        match c.__name__:
            case "ContinuousGalerkinMixin":    return cls.ContinuousGalerkin
            case "DiscontinuousGalerkinMixin": return cls.DiscontinuousGalerkin
            case _: return None

    def __str__(self) -> str:
        match self:
            case self.ContinuousGalerkin:    return "continous galerkin"
            case self.DiscontinuousGalerkin: return "discontinuous galerkin"

    def __repr__(self) -> str: return str(self)

    def abbrev(self) -> str:
        match self:
            case self.ContinuousGalerkin:    return "cg"
            case self.DiscontinuousGalerkin: return "dg"
        

def tests():
    scheme = DiscretisationScheme.ContinuousGalerkin
    print(scheme)

if __name__ == "__main__":
    tests()