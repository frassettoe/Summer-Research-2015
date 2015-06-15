#Michael and Erin
#06/15/15
#Test of optimization function, using only function calls

import numpy
from scipy.optimize import minimize

def simpleFunction(a,b,c):
    d = 5*a**2+3*b**2+c**2
    return d


class testDevice:
    def simplerFunction(self,a):
        b = self.part1(a)
        b = self.part2(b)
        return b

    def part1(self,a):
        return a**2
    def part2(self,a):
        return 3*a

    def opt(self):
        res = minimize(self.simplerFunction, self.x0 ,method = 'Newton-CG',jac = simplerFunctionDer,hessp = simplerFunctionHes,options={'disp':True})
        print(res.x)
        return res

    def __init__(self):
        self.x0 = numpy.array([1000])


def simplerFunctionDer(a):
    return 6*a

def simplerFunctionHes(a,b):
    return 6

obj = testDevice()
print(obj.opt())