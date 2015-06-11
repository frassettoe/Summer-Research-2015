__author__ = 'Owner'

import math
import random
import copy
import numpy
from scipy.optimize import minimize
import scipy

#input: coformal variations
#output: edgelength of edge ij
def EdgeLength(i,j):
    edgelength= math.exp(.5*(i+j))
    return edgelength

#input: conformal variations
#output: the face angle i in triangle i,j,k
def FaceAngleijk(i, j, k): #i,j,k
    angle = math.acos((EdgeLength(i,j)**2+EdgeLength(i,k)**2-EdgeLength(j,k)**2)/(2*EdgeLength(i,j)*EdgeLength(i,k)))
    return angle

#input:conformal variations
#outpu: the dihedral angle around edge ij in tetrahedra ijkl
def DihedralAngleijkl(i,j,k,l):
    diangle = math.acos((math.cos(FaceAngleijk(i,k,l))-math.cos(FaceAngleijk(i,j,k))*math.cos(FaceAngleijk(i,j,l)))/(math.sin(FaceAngleijk(i,j,k))*math.sin(FaceAngleijk(i,j,l))))
    return diangle


#input: conformal variations
#output: edge curvature of edge ij, in tetrahedra ijkl
def EdgeCuvature(i,j,k,l):
    ecurvature = (2*math.pi-DihedralAngleijkl(i,j,k,l))*EdgeLength(i,j)
    return ecurvature

#input: conformal variations
#output: curvature around vertex i
def VertexCurvature(i,j,k,l):
    vcurvature=.5*(EdgeCuvature(i,j,k,l)+EdgeCuvature(i,k,j,l)+EdgeCuvature(i,l,j,k))
    return vcurvature

#input:conformal variation
#output: the sum of all of the edges
def TotalEdgeLength(i,j,k,l):
    totallength = EdgeLength(i,j)+EdgeLength(i,k)+EdgeLength(i,l)+EdgeLength(k,j)+EdgeLength(k,l)+EdgeLength(j,l)
    return totallength

#input: conformal variation
#output: LEHR
def LEHR(i,j,k,l):
    LEHR= (VertexCurvature(i,j,k,l)+VertexCurvature(j,i,k,l)+VertexCurvature(k,i,j,l)+VertexCurvature(l,i,j,k))/TotalEdgeLength(i,j,k,l)
    return LEHR

#input:
#output:
def gradient(i,j,k,l):
    Grad =
    return Grad

print(EdgeLength(.5,.5))
print(FaceAngleijk(.5,.5,.5))
print(DihedralAngleijkl(.5,.5,.5,.5))
print(EdgeCuvature(.5,.5,.5,.5))
print(VertexCurvature(.4,.5,.3,.5))
print(VertexCurvature(.5,.4,.3,.5))
print(TotalEdgeLength(.5,.5,.5,.5))
print(LEHR(.5,.5,.5,.5))