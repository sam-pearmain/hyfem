import abc
import ufl

from hyfem.core.spaces import Spaces
from hyfem.core.eqns.traits import LinearSolveable, NonlinearSolveable
from hyfem.utils import *

class LinearEquation(LinearSolveable):
    """
    A closed-form equation with a single unknown variable
    """
    def __init__(self) -> None:
        super().__init__()

class NonlinearEquation(NonlinearSolveable):
    ...