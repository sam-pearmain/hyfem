__all__ = [
    "TodoError", 
    "todo"
]

class TodoError(NotImplementedError): ...

def todo(*args) -> None: raise TodoError(*args)


def tests():
    todo("we need to do this") 

if __name__ == "__main__":
    tests()