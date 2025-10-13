import functools
from typing import Any, Callable, TypeVar

__any__ = [
    "Property", 
    "ClassMethod", 
    "StaticMethod",
    "error", 
    "check", 
]

Property = property
ClassMethod = classmethod
StaticMethod = staticmethod

E = TypeVar('E', bound = Exception)
def error(func: Callable[..., E]) -> Callable[..., E]:
    @functools.wraps(func)
    def wrapper(*args: Any, **kwds: Any) -> E:
        result = func(*args, **kwds)
        if not isinstance(result, Exception):
            raise TypeError(
                f"{func.__name__} must return an Exception but returned {type(result).__name__}"
            )
        return result
    return wrapper

def check(func: Callable[..., bool]) -> Callable[..., bool]:
    @functools.wraps(func)
    def wrapper(*args: Any, **kwds: Any) -> bool:
        result = func(*args, **kwds)
        if not isinstance(result, bool):
            raise TypeError(
                f"{func.__name__} must return a bool but returned {type(result).__name__}"
            )
        return result
    return wrapper