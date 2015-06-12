__author__ = 'Owner'

import math
import random
import copy
import numpy
from scipy.optimize import minimize
import scipy


#inport global variable background lengths
#accesss to triangle list and tetrahedralist
def DefineGlobalVariables(manifoldFile = 'manifoldExample.txt'):
    readFile = open(manifoldFile)
    data = readFile.read()        #Prepares file for read in
    data = data.split("facets :=") #Look up strip to remove white space
    data[1] = data[1].strip('[];')
    data[1] = data[1].split('],[')
    tetrahedra = []
    for i in range(0, len(data[1])):   #List comprehensions
        tetrahedra.append(data[1][i])
    for i in range(0, len(tetrahedra)):
        tetrahedra[i] = tetrahedra[i].split(',')
    readFile.close()
    tetrahedra = [[int(i) for i in tetrahedra[j]] for j in range(len(tetrahedra))]

#input: coformal variations
#output: edgelength of edge ij
def EdgeLength(i,j):
    edgelength= math.exp(.5*(i+j))#times the background length
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
    ecurvature = (2*math.pi-DihedralAngleijkl(i,j,k,l))*EdgeLength(i,j) # missing a sum of the dihedral angles, should be the sum of the dihedral angle of that edge in all tetrahedra that edge appears in
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

#inptut: confromal variations
#output: hij,k which is used in the lij* calculations in the Hessian
def hijk(i,j,k):
    hijk= (EdgeLength(i,k)-EdgeLength(i,j)*math.cos(FaceAngleijk(i,j,k)))/(math.sin(FaceAngleijk(i,j,k)))
    return hijk

#input: conformal variations
#output: hijk,l which is used in the lij* calculations in the Hessian
def hijkl(i,j,k,l):
    hijkl=(hijk(i,j,l)-hijk(i,j,l)*math.cos(DihedralAngleijkl(i,j,k,l)))/(math.sin(DihedralAngleijkl(i,j,k,l)))
    return hijkl

#input: conformal variations
#output: lij* which is used in the Hessian
def lijstar(i,j,k,l):
    lijstar=.5*(hijk(i,j,k)*hijkl(i,j,k,l)+hijk(i,j,l)*hijkl(i,j,l,k)) # missing a sum of all of the tetrahedra the edge ij is in which goes infromt
    return lijstar

#input:
#output:
def gradient(i,j,k,l):
    Grad = Grad
    return Grad

print(EdgeLength(.5,.5))
print(FaceAngleijk(.5,.5,.5))
print(DihedralAngleijkl(.5,.5,.5,.5))
print(EdgeCuvature(.5,.5,.5,.5))
print(VertexCurvature(.4,.5,.3,.5))
print(VertexCurvature(.5,.4,.3,.5))
print(TotalEdgeLength(.5,.5,.5,.5))
print(LEHR(.5,.5,.5,.5))