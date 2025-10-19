from test import Equation


class Burgers(Nonlinear, TimeDependent, PDE):
    def state_variable(self) -> str:
        return ['u']