import functools
from typing import Any, Callable, TypeVar

__any__ = [
    "String", 
    "Int", 
    "Float", 
    "Property", 
    "ClassMethod", 
    "StaticMethod",
]

# types
String = str
Int = int
Float = float

# decors
Property = property
ClassMethod = classmethod
StaticMethod = staticmethod