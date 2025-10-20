import abc

from hyfem.core.spaces import Spaces

class Solvable(abc.ABC):
    """The root base class for any equation or system of equations"""
    _spaces: Spaces | None

    def __init__(self) -> None:
        super().__init__()
        self._spaces = None

    def assign_function_spaces(self, spaces: Spaces) -> None:
        supported_eqn, _ = spaces.defined_on_str()
        if not supported_eqn == type(self).__name__:
            raise ValueError(f"attempted to assign <Spaces ({supported_eqn})> object to {type(self).__name__}")
        self._spaces = spaces

    def has_function_spaces(self) -> bool:
        return self._spaces is not None

    def has_defined_function_spaces(self) -> bool:
        if not self.is_finalised():
            raise RuntimeError(f"{type(self).__name__} has not been assigned a {Spaces.__name__} object")
        return self._spaces.all_spaces_assigned()

class Unclosed(abc.ABC):
    """A marker trait for unclosed PDEs"""
    ...

class Spatial(abc.ABC):
    """A trait for PDEs with a spatial dimension"""
    ...

class Temporal(abc.ABC):
    """A trait for PDEs with a temporal dimension"""
    ...

class LinearInSpace(Spatial):
    """A trait for PDEs that are linear in its spatial domain"""
    ...

class NonlinearInSpace(Spatial):
    """A trait for PDEs that are nonlinear in its spatial domain"""
    ...

class LinearInTime(Temporal):
    """A trait for PDEs which are linear in time"""
    ...

class NonlinearInTime(Temporal):
    """A trait for PDEs which are nonlinear in time"""
    ...
