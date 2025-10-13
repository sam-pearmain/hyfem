class TodoError(NotImplementedError): ...

def todo(*args) -> None:
    """Raise a TodoError"""
    raise TodoError(*args)


def tests():
    todo("we need to do this") 

if __name__ == "__main__":
    tests()