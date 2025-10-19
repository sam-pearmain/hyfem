import abc

from hyfem.utils import *

class PDE(abc.ABC):
    """
    A closed-form PDE with a single unknown variable
    """

    @Property
    def unknown(self) -> str: return self._unknown_impl()

    def is_system(self) -> bool: return False

    abc.abstractmethod
    def _unknown_impl(self) -> str: ...