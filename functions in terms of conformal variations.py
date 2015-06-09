__author__ = 'Owner'

import math
import random
import copy
import numpy
from scipy.optimize import minimize
import scipy

#input: conformal variations
#output: the face angle of vertex with convar1 opposite edge with conformal variations 2,3
def FaceAngle(convar1, convar2, convar3): #i,j,k
    angle = math.acos((math.exp(.5*(convar1+convar2))**2+math.exp(.5*(convar1+convar3))**2-math.exp(.5*(convar2+convar3))**2)/(2*math.exp(.5*(convar1+convar2))*math.exp(.5*(convar1+convar3))))
    return angle

#input:conformal variations
#outpu: the dihedral angle
def DihedralAngle(convar1,convar2,convar3,convar4):
    diangle = math.acos((math.cos(FaceAngle(convar1,convar3,convar4))-math.cos(FaceAngle(convar1,convar2,convar3))*math.cos(FaceAngle(convar1,convar2,convar4)))/(math.sin((FaceAngle(convar1,convar2,convar3))*math.sin(FaceAngle(convar1,convar2,convar4)))))
    print(math.cos(FaceAngle(convar1,convar3,convar4)))
    print(math.cos(FaceAngle(convar1,convar2,convar3)))
    print(math.cos(FaceAngle(convar1,convar2,convar4)))
    print(math.sin(FaceAngle(convar1,convar2,convar3)))
    print(math.sin(FaceAngle(convar1,convar2,convar4)))
    print(diangle)
    return diangle

print(FaceAngle(.5,.5,.5))
print(DihedralAngle(.5,.5,.5,.5))