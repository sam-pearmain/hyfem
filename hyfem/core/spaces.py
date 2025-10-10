from typing import List, Tuple


class Spaces(object):
    _spaces: List[Tuple[str, int]]
    
    def __init__(self) -> None:
        self._spaces = []

    def add_function_space(self, family: str, degree: int) -> None:
        if (family, degree) not in self._spaces:
            self._spaces.append((family, degree))