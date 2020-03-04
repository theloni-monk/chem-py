import math

#TODO: WRITEME sciNum
class sciNum:    
    def __init__(self, base, exponent):
        self.base = base
        self.exp = exponent
        self.sigfigs = _getSigFigs(base)

    @classmethod
    def fromString(cls, str):
        pass

    @classmethod
    def fromLong(cls, long):
        pass

    def __str__(self):
        pass

    def __add__(self, other):
        pass

    def __sub__(self, other):
        pass

    def __mul__(self, other):
        pass

    def __div__(self, other):
        pass

def _getSigFigs(num):
    pass