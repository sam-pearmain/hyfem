from hyfem.equations.base import Equation


class Advection(Equation):
    def state_variables(cls):
        return super().state_variables()