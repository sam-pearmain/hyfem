import abc
from typing import List

class Equation(abc.ABC):
    @classmethod
    @abc.abstractmethod
    def state_variables(cls) -> List[str]: 
        """The state variables of the equation, for example [u, rho_u, E]"""
        ...

    @classmethod
    @abc.abstractmethod
    def auxiliary_variables(cls) -> List[str]: 
        """The auxiliary variables of the equation, for example [k, ε]"""
        ...

    @classmethod
    def n_state_variables(cls) -> int: return len(cls.state_variables())

    @classmethod
    def n_auxiliary_variables(cls) -> int: return len(cls.auxiliary_variables())