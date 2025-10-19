from typing import List, Type

from hyfem.utils import *


class Spaces:
    _vars: List[str]

    def __init__(self, eqn: Type[PDE] | Type[System]):
        if isinstance(eqn):
            raise 
        
        self._vars = eqn.state_variables()
