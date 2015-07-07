__author__ = 'Owner'

# uncomment all lines with 333 in them to get hessian working

#tp://www.ripon.edu/wp-content/uploads/2012/10/doubletetra.pdf
#10/16/14 Michael: Added advancedCheckLegalTetrahedra function for checking legal tetrahedron and used it main, made minor change checkLegalTetrahedron
import math
import random
import copy
import numpy
from scipy.optimize import minimize
import scipy
import time
import cProfile
import BasisBuilder




class Edge:
    #vertex1: an integer that corresponds with the first part of the name of an edge
    #vertex2: an integer that corresponds with the second part of the name of an edge
    #edgelength: a number that represents how long the edge is
    #edgecurvature: a number that is edge curvature
    #tetrahdraEdgeIsIn: location in tetrahedraList of all tetrahedra that contain that edge

    #__init__(self,vertex1 = 1,vertex2 = 2): returns nothing, sets default values
    #calculateEdgeCurvature(self,listOfTetrahedra): returns nothing, sets the edge curvature

    #input: three numbers that represent the name of an edge and the edge length respectively
    #output: none, sets default values of members of the edge class
    #Author: Prof. Young, 10/2/2014
    #change log: Michael  10/27/2014
    def __init__(self,vertex1 = 1,vertex2 = 2, length = 1):
        self.vertex1 = vertex1
        self.vertex2 = vertex2
        self.edgelength = length
        self.edgecurvature = 0
        self.tetrahedraEdgeIsIn = [] #This a lis of numbers that corresponds to the location of the tetrahedron in a list of Tetrahedron (attribute of metric class and backgroundMetricClass class)
        self.edgelengthStar = 0  #Used in Hessian calculation
    #10/3/2014: added tetrahedraEdgeIsIn, MELT
    #10/27/14: Vertex class now takes a third argument, length, with default value of one


    #input: list of Tetrahedron objects
    #output: none, sets edgecurvature to the curvature of the edge
    #author: MELT, 10/14/2014
    #change log: none
    def calculateEdgeCurvature(self,listOfTetrahedra):
        diList = []  #empty list of dihedral angles
        for i in range (len(self.tetrahedraEdgeIsIn)):  #for each tetrahedran the edge is in
            tetLocation = self.tetrahedraEdgeIsIn[i]    #find the the location of that tetrahedron in the list of tetrahedra
            # Finds the dihedral angle of edge for that tetrahedra
            #Dihedral angle N at position L corresponds with position of Edge N at position L in edgesintetrahedron
            singleDi = listOfTetrahedra[tetLocation].dihedralanglelist[listOfTetrahedra[tetLocation].edgesintetrahedron.index([self.vertex1,self.vertex2])]
            diList.append(singleDi)  #ad dihedral angle to list of dihedral angles
        self.edgecurvature = (2*math.pi-sum(diList))*self.edgelength # Dihedral angle formula
    #


    #input: list of Tetrahedron objects
    #output: none, sets edgelengthStar to the correct edge length
    #author: ME, 02/02/2015
    #change log: none
    def calculateEdgeLengthStar(self, listOfTetrahedra):
        starList = []
        for i in range (len(self.tetrahedraEdgeIsIn)):
            tetLocation = self.tetrahedraEdgeIsIn[i]
            list  = copy.deepcopy(listOfTetrahedra[tetLocation].vertexList)
            list.sort()
            iLoc = list[list.index(self.vertex1)]
            jLoc = list[list.index(self.vertex2)]
            list.remove(self.vertex1)
            list.remove(self.vertex2)
            kLoc = list[0]
            lLoc = list[1]
            lface = listOfTetrahedra[tetLocation].vertexList.index(lLoc)
            kface = listOfTetrahedra[tetLocation].vertexList.index(kLoc)
            hijk = listOfTetrahedra[tetLocation].faceList[lface].hList[listOfTetrahedra[tetLocation].faceList[lface].vertexList.index(kLoc)]
            hijkl = listOfTetrahedra[tetLocation].tetCenDisList[lface]
            hijl = listOfTetrahedra[tetLocation].faceList[kface].hList[listOfTetrahedra[tetLocation].faceList[kface].vertexList.index(lLoc)]
            hijlk = listOfTetrahedra[tetLocation].tetCenDisList[kface]
            singleStar = hijk*hijkl+hijl*hijlk
            starList.append(singleStar)
        self.edgelengthStar = .5*sum(starList)




class face:
    def getAngle(self,c,a,b):  #uses law of cosines to find the angle given three sides of a triangle, returns cosine of angle
        temp = a**2+b**2-c**2
        temp = temp/(2*a*b)
        return temp
    #
    def getTriCenDis(self,leg2,leg1,Angle):
        temp=(.5*leg2-.5*leg1*math.cos(Angle))/math.sin(Angle)
        return temp

    def __init__(self, edgeTable, vertex1 = 1, vertex2 = 2, vertex3 = 3):
        self.vertexList = [vertex1,vertex2,vertex3]
        self.edgeLength = [edgeTable[vertex2][vertex3].edgelength,edgeTable[vertex1][vertex3].edgelength,edgeTable[vertex1][vertex2].edgelength]  #edgeLength[i] is opposite to vertex i-1
        edge1 = self.edgeLength[0]
        edge2 = self.edgeLength[1]
        edge3 = self.edgeLength[2]
        self.angleList = [math.acos(self.getAngle(edge1,edge2,edge3)),math.acos(self.getAngle(edge2,edge1,edge3)),math.acos(self.getAngle(edge3,edge2,edge1))]
        #self.hList = [self.getTriCenDis(edge3,edge1,self.angleList[1]),self.getTriCenDis(edge1,edge2,self.angleList[2]),self.getTriCenDis(edge2,edge3,self.angleList[0])]   #333 used in hessian


class Tetrahedron:
    #vertex1: an integer that corresponds to the first part of the name of a tetrahedron
    #vertex2: an integer that corresponds to the second part of the name of a tetrahedron
    #vertex3: an integer that corresponds to the third part of the name of a tetrahedron
    #vertex4: an integer that corresponds to the fourth part of the name of a tetrahedron
    #edgesintetrahedron: list of lists of names of edges, edge 1,2 is stored as [1,2]
    #dihedralanglelist: list of dihedral angles where dihedralanglelist[i] is the dihedral angle of the edge with name at edgesintetrahedron[i]

    #__init__(self,vertex1 = 1, vertex2 = 2, vertex3 = 3, vertex4 = 4): returns nothing, sets default values
    #checkLegalTetrahedron(self,tableOfEdges): returns true or false, if true all triangle inequalities and CM determinants are legit, if false not legit
    #getAngle(self,c,a,b): returns temp, the cosine of the angle opposite of edge c; mini-function for calculateDihedralAngles
    #calDiAngle(self,Eij,Eik,Eil,Ejk,Ejl,Ekl): returns temp, applies dihedral angle formula, mini-function for calculateDihedralAngles
    #calculateDihedralAngles(self,tableOfEdges): returns nothing, mega-function that finds the dihedral angles of every edge and puts them in the dihedralanglelist

    #input: four numbers that represent the name of a tetrahedron
    #output: none, sets default values of members of the tetrahedron class
    #author:Prof. Young, 10/2/14
    #change log: ?
    def __init__(self,vertex1 = 1, vertex2 = 2, vertex3 = 3, vertex4 = 4):
        self.vertex1 = vertex1
        self.vertex2 = vertex2
        self.vertex3 = vertex3
        self.vertex4 = vertex4
        self.vertexList = [vertex1,vertex2,vertex3,vertex4]
        self.edgesintetrahedron = []  #edges in the tetrahedron stored as lists ex, edge 1,2 stored as [1,2]
        self.dihedralanglelist = []
        self.tetCenDisList = [0]*4
        self.faceList = [0]*4 #faceList[i] is the face oppositve vertex i-1
        self.angleTable = []
    #


    def initalizeTriangles(self,tableOfEdges):
        self.faceList[0] = face(tableOfEdges,self.vertex2,self.vertex3,self.vertex4)
        self.faceList[1] = face(tableOfEdges,self.vertex1,self.vertex3,self.vertex4)
        self.faceList[2] = face(tableOfEdges,self.vertex1,self.vertex2,self.vertex4)
        self.faceList[3] = face(tableOfEdges,self.vertex1,self.vertex2,self.vertex3)
    def checkTriangles(self,edge1,edge2,edge3,edge4,edge5,edge6):
        legalTriangle = True
        if edge1+edge2+edge4-2*max(edge1,edge2,edge4) <= 0:
            legalTriangle=False
        if edge1+edge3+edge5-2*max(edge1,edge3,edge5) <= 0:
            legalTriangle=False
        if edge2+edge3+edge6-2*max(edge2,edge3,edge6) <= 0:
            legalTriangle=False
        if edge4+edge5+edge6-2*max(edge4,edge5,edge6) <= 0:
            legalTriangle=False
        return legalTriangle
    checkLCalls = 0
    TetBreak = 0
    TriBreak = 0
    BothBreak = 0
    #input: Table of edge objects
    #output: finalDecision, a boolean that tells whether legality conditions have been met
    #author: MELT, 10/13/2014
    #change log: Michael 10/16/14
    def checkLegalTetrahedron(self,tableOfEdges):
        Tetrahedron.checkLCalls = Tetrahedron.checkLCalls+1
        legalTriangle=True
        legalTetrahedron=True
        edge1=tableOfEdges[self.vertex1][self.vertex2].edgelength
        edge2=tableOfEdges[self.vertex1][self.vertex3].edgelength
        edge3=tableOfEdges[self.vertex1][self.vertex4].edgelength
        edge4=tableOfEdges[self.vertex2][self.vertex3].edgelength
        edge5=tableOfEdges[self.vertex2][self.vertex4].edgelength
        edge6=tableOfEdges[self.vertex3][self.vertex4].edgelength
        det = (-2)*(edge1**4)*(edge6**2)-(2)*(edge1**2)*(edge2**2)*(edge4**2)+(2)*(edge1**2)*(edge2**2)*(edge5**2)+(2)*(edge1**2)*(edge2**2)*(edge6**2)+(2)*(edge1**2)*(edge3**2)*(edge4**2)-(2)*(edge1**2)*(edge3**2)*(edge5**2)+(2)*(edge1**2)*(edge3**2)*(edge6**2)+(2)*(edge1**2)*(edge4**2)*(edge6**2)+(2)*(edge1**2)*(edge5**2)*(edge6**2)-(2)*(edge1**2)*(edge6**4)-(2)*(edge2**4)*(edge5**2)+(2)*(edge2**2)*(edge3**2)*(edge4**2)+(2)*(edge2**2)*(edge3**2)*(edge5**2)-(2)*(edge2**2)*(edge3**2)*(edge6**2)+(2)*(edge2**2)*(edge4**2)*(edge5**2)-(2)*(edge2**2)*(edge5**4)+(2)*(edge2**2)*(edge5**2)*(edge6**2)-(2)*(edge3**4)*(edge4**2)-(2)*(edge3**2)*(edge4**4)+(2)*(edge3**2)*(edge4**2)*(edge5**2)+(2)*(edge3**2)*(edge4**2)*(edge6**2)-(2)*(edge4**2)*(edge5**2)*(edge6**2)
        if det<=0:
            legalTetrahedron=False
            return False
        legalTriangle = self.checkTriangles(edge1,edge2,edge3,edge4,edge5,edge6)
        # if legalTetrahedron == False and legalTriangle == False:
        #     Tetrahedron.BothBreak = Tetrahedron.BothBreak+1
        # elif legalTetrahedron == False:
        #     Tetrahedron.TetBreak = Tetrahedron.TetBreak+1
        # elif legalTriangle == False:
        #     Tetrahedron.TriBreak = Tetrahedron.TriBreak+1
        # print("BOT: "+str(Tetrahedron.BothBreak))
        # print("TET: "+str(Tetrahedron.TetBreak))
        # print("DEL: "+str(Tetrahedron.TriBreak))
        finalDecision=legalTetrahedron and legalTriangle
        return finalDecision
    # 10/16/14; Michael; change return(finalDecision) to return finalDecision

    #input: Table of edge objects
    #output: volume of the terahedron
    #author: Erin 6/29/2015
    #change log:
    def volumeOfTetrahedron(self,tableOfEdges):
        edge1=tableOfEdges[self.vertex1][self.vertex2].edgelength
        edge2=tableOfEdges[self.vertex1][self.vertex3].edgelength
        edge3=tableOfEdges[self.vertex1][self.vertex4].edgelength
        edge4=tableOfEdges[self.vertex2][self.vertex3].edgelength
        edge5=tableOfEdges[self.vertex2][self.vertex4].edgelength
        edge6=tableOfEdges[self.vertex3][self.vertex4].edgelength
        det = (-2)*(edge1**4)*(edge6**2)-(2)*(edge1**2)*(edge2**2)*(edge4**2)+(2)*(edge1**2)*(edge2**2)*(edge5**2)+(2)*(edge1**2)*(edge2**2)*(edge6**2)+(2)*(edge1**2)*(edge3**2)*(edge4**2)-(2)*(edge1**2)*(edge3**2)*(edge5**2)+(2)*(edge1**2)*(edge3**2)*(edge6**2)+(2)*(edge1**2)*(edge4**2)*(edge6**2)+(2)*(edge1**2)*(edge5**2)*(edge6**2)-(2)*(edge1**2)*(edge6**4)-(2)*(edge2**4)*(edge5**2)+(2)*(edge2**2)*(edge3**2)*(edge4**2)+(2)*(edge2**2)*(edge3**2)*(edge5**2)-(2)*(edge2**2)*(edge3**2)*(edge6**2)+(2)*(edge2**2)*(edge4**2)*(edge5**2)-(2)*(edge2**2)*(edge5**4)+(2)*(edge2**2)*(edge5**2)*(edge6**2)-(2)*(edge3**4)*(edge4**2)-(2)*(edge3**2)*(edge4**4)+(2)*(edge3**2)*(edge4**2)*(edge5**2)+(2)*(edge3**2)*(edge4**2)*(edge6**2)-(2)*(edge4**2)*(edge5**2)*(edge6**2)
        if det>0:
            volume = math.sqrt(det/288)
        else:
            det=0
            volume = math.sqrt(det/288)
        return volume


    #input: c,a,b, which are edge lengths
    #output: temp, the cosine of the angle opposite of edge c
    #author: MELT, 10/9/2014
    #change log:
    def getAngle(self,c,a,b):  #Uses law of cosines to return the cosine of an angle given three edge lengths of a triangle
        temp = a**2+b**2-c**2
        temp = temp/(2*a*b)
        return temp
    #


    #input: leg1A and leg2A are connected to target vertex and Target is between these two legs.  leg1B and leg2B are remaining lets not connected to target vertex.  Dihedral angle is the dihedral anlge of the unused leg
    #output: TetCenDis and related values
    #author: ME, 02/02/2015
    #change log:
    def calTetCenDis(self,i,j,k,l):
         list = [i,j,l]
         list.sort()
         for m in range(len(self.faceList)):
             if self.faceList[m].vertexList == list:
                for e in range(len(self.faceList[m].hList)):
                  if self.faceList[m].vertexList[e] == l:
                      hijl = self.faceList[m].hList[e]
                      break
         list = [i,j,k]
         list.sort()
         for m in range(len(self.faceList)):
             if self.faceList[m].vertexList == list:
                for e in range(len(self.faceList[m].hList)):
                    if self.faceList[m].vertexList[e] == k:
                        hijk = self.faceList[m].hList[e]
                        break
         list = [i,j]
         list.sort()
         DihedralAngleLocation = self.edgesintetrahedron.index(list)
         DihedralAngleijkl = self.dihedralanglelist[DihedralAngleLocation]
         result = (hijl-hijk*math.cos(DihedralAngleijkl))/math.sin(DihedralAngleijkl)
         return result

    def getTetCenDis(self):
        self.tetCenDisList[0] = self.calTetCenDis(self.vertex2,self.vertex3,self.vertex4,self.vertex1)
        self.tetCenDisList[1] = self.calTetCenDis(self.vertex1,self.vertex3,self.vertex4,self.vertex2)
        self.tetCenDisList[2] = self.calTetCenDis(self.vertex2,self.vertex4,self.vertex1,self.vertex3)
        self.tetCenDisList[3] = self.calTetCenDis(self.vertex3,self.vertex2,self.vertex1,self.vertex4)

    #input: Eij,Eik,Eil,Ejk,Ejl,Ekl, numbers that represent edge lengths of a tetrahedron, edge Eij is the dihedral angle
    #output: temp, result of applying dihedral angle formula
    #author: MELT, 10/9/2014
    #change log:
    def calDiAngle(self,Eij,Eik,Eil,Ejk,Ejl,Ekl):  #Used to find a single dihedral angle
        temp = self.getAngle(Ekl,Eik,Eil)-self.getAngle(Ejk,Eij,Eik)*self.getAngle(Ejl,Eij,Eil)
        temp = temp/math.sin(math.acos(self.getAngle(Ejk,Eij,Eik)))
        temp = temp/math.sin(math.acos(self.getAngle(Ejl,Eij,Eil)))
        return temp
    #Should find dihedral angle of edge ij
    def calDiAngle2(self,Ailk,Aijk,Aijl):  #Used to find a single dihedral angle
        temp = Ailk-Aijk*Aijl
        temp = temp/math.sin(math.acos(Aijk))
        temp = temp/math.sin(math.acos(Aijl))
        return temp

        #temp = angle oposite Edge kl-angle oposite edge kj * angle opposite edge jl
    #

    #input: table of edge objects
    #output: none, finds the dihedral angle at each edge and adds it to the dihedralanglelist
    #author: MELT, 10/9/2014
    #change log:
    def calculateDihedralAngles(self,tableOfEdges):  #Used to find all the dihedral angles of a tetrahedron
        self.dihedralanglelist = []
        angleTable = [0]*12
        edge12=tableOfEdges[self.vertex1][self.vertex2].edgelength  #gets edge lengths in a form easier for the coder to understand
        edge13=tableOfEdges[self.vertex1][self.vertex3].edgelength  #this reduced the risk of mistyping the order of edges in calDiAngle calls
        edge14=tableOfEdges[self.vertex1][self.vertex4].edgelength
        edge23=tableOfEdges[self.vertex2][self.vertex3].edgelength
        edge24=tableOfEdges[self.vertex2][self.vertex4].edgelength
        edge34=tableOfEdges[self.vertex3][self.vertex4].edgelength
        angleTable[0] = self.getAngle(edge12,edge23,edge13)
        angleTable[1] = self.getAngle(edge23,edge12,edge13)
        angleTable[2] = self.getAngle(edge13,edge23,edge12)
        angleTable[3] = 0#self.getAngle(edge13,edge14,edge34)
        angleTable[4] = self.getAngle(edge14,edge13,edge34)
        angleTable[5] = self.getAngle(edge34,edge14,edge13)
        angleTable[6] = 0#self.getAngle(edge12,edge14,edge24)
        angleTable[7] = self.getAngle(edge14,edge12,edge24)
        angleTable[8] = self.getAngle(edge24,edge14,edge12)
        angleTable[9] = self.getAngle(edge24,edge34,edge23)
        angleTable[10] = self.getAngle(edge34,edge24,edge23)
        angleTable[11] = -1#self.getAngle(edge23,edge34,edge24)
        self.dihedralanglelist.append(math.acos(self.calDiAngle2(angleTable[5],angleTable[1],angleTable[8])))  #finds each dihedral angle
        self.dihedralanglelist.append(math.acos(self.calDiAngle2(angleTable[8],angleTable[5],angleTable[1])))
        self.dihedralanglelist.append(math.acos(self.calDiAngle2(angleTable[1],angleTable[8],angleTable[5])))
        self.dihedralanglelist.append(math.acos(self.calDiAngle2(angleTable[7],angleTable[2],angleTable[10])))
        self.dihedralanglelist.append(math.acos(self.calDiAngle2(angleTable[2],angleTable[10],angleTable[7])))
        self.dihedralanglelist.append(math.acos(self.calDiAngle2(angleTable[0],angleTable[9],angleTable[4])))
        # print(self.edgesintetrahedron[0],self.dihedralanglelist[0],self.vertex1,self.vertex2)  #prints dihedral angles, used only for checking for errors
        # print(self.edgesintetrahedron[1],self.dihedralanglelist[1],self.vertex1,self.vertex3)
        # print(self.edgesintetrahedron[2],self.dihedralanglelist[2],self.vertex1,self.vertex4)
        # print(self.edgesintetrahedron[3],self.dihedralanglelist[3],self.vertex2,self.vertex3)
        # print(self.edgesintetrahedron[4],self.dihedralanglelist[4],self.vertex2,self.vertex4)
        # print(self.edgesintetrahedron[5],self.dihedralanglelist[5],self.vertex3,self.vertex4)
        # for i in range(len(self.dihedralanglelist)):
        #     print(self.dihedralanglelist[i])
    #



class backgroundMetricClass: #need to change so conformal variations are always 0.
    #input: list of edge names and table of edge objects
    #output: a number that is LEHR
    #author: MELT, 10/9/2014
    #change log:
    def findLEHR(self,listOfEdges,tableOfEdges):
        listOfLengths = []
        listOfCurvatures = []
        for i in range(len(listOfEdges)):  #for each edge, adds length and curvature of edge to respective list
            edgeSpot1 = listOfEdges[i][0]
            edgeSpot2 = listOfEdges[i][1]
            listOfLengths.append(tableOfEdges[edgeSpot1][edgeSpot2].edgelength)
            listOfCurvatures.append(tableOfEdges[edgeSpot1][edgeSpot2].edgecurvature)
        return sum(listOfCurvatures)/sum(listOfLengths)  #preforms LEHR formula
    #

    #input: empty table of edges, number of vertices
    #output: none, makes the edge table of proper size, i.e. 15+1
    #author: MELT, 10/3/14
    #change log: M 10/3/14
    def createEdgeTable(self,edgetable,numberOfVertices = 15):
        numberOfVertices = numberOfVertices+1  #adds one, as this way edge 7 is at position 7, not position 6
        for row in range(numberOfVertices):    #(nessessary due to lists starting at position 0)
            edgetable.append([])
            for column in range(numberOfVertices):
                edgetable[row].append(0)  #sets every edge location to 0, indicating that no edges are present
    #10/3/14 M, edgetable is now an input

    #input: list of tetrahedron object, table of edge objects, list of edge names, file of background metric created from LEHRBackgroundMetricBuilder, word NONE assumes all edge lengths one and no file is real
    #output: none, fills the edge table with edge objects and puts the tetrahedron into the list of tetrahedraEdgeIsIn
    #author: MELT, 10/6/2014
    #change log: 11/18/14 ME
    def fillEdgeTable(self,listOfTetrahedra,tableOfEdges,listOfEdges,fileName = "backgroundMetric.txt"):
        #Reads in a background metric file and forms the file into a useful format
        backgroundMetric = open(fileName,"r")
        storage = backgroundMetric.readlines()
        nameList = []
        lengthList = []
        i = 0
        while(i < len(storage)-1):
            #breaks background metric into a names and lengths, such that name at position i corresponds to length and position i
            nameList.append(storage[i])
            lengthList.append(float(storage[i+1]))
            i = i+2
        # makes files into useful format
        for i in range(len(nameList)):
            nameList[i] = nameList[i].split(",")
            nameList[i][0] = int(nameList[i][0])
            nameList[i][1] = int(nameList[i][1])
        counter = 0
        for i in range(len(listOfTetrahedra)):  #for each tetrahedra
            for j in range(len(listOfTetrahedra[i].edgesintetrahedron)):  #look at every edge in the tetrahedra
                    # If edge is not in the the tableOfEdges, add it
                    edgeInTetrahedronIndex1 = listOfTetrahedra[i].edgesintetrahedron[j][0]
                    edgesInTetrahedronIndex2 = listOfTetrahedra[i].edgesintetrahedron[j][1]
                    if tableOfEdges[edgeInTetrahedronIndex1][edgesInTetrahedronIndex2] == 0:
                        # Assigns edge its name
                        tableOfEdges[edgeInTetrahedronIndex1][edgesInTetrahedronIndex2] = Edge(edgeInTetrahedronIndex1,edgesInTetrahedronIndex2,lengthList[counter])
                        counter = counter+1  #increases the number of edges counted by one
                        # Adds tetrahedran to list of tetrahedra edge is in
                        tableOfEdges[edgeInTetrahedronIndex1][edgesInTetrahedronIndex2].tetrahedraEdgeIsIn.append(i)
                        listOfEdges.append([edgeInTetrahedronIndex1,edgesInTetrahedronIndex2])
                    # adds edge to edgesInTetrahedra list if edge already exists in the table
                    else:
                        tableOfEdges[edgeInTetrahedronIndex1][edgesInTetrahedronIndex2].tetrahedraEdgeIsIn.append(i)
    #10/29/14: Michael Added new ability! Function can now be given a backgroundmetric and a list of conformal variations.
    #Both must be given in order for new function attibutes to be used.  Old use without conformal variations and background metric can still be used.
    #11/18/14: Micahel and Erin removed the ability to not take a background metric and conformal variations.

    #Input: The two vertex conformal varations and the length of the edge in question
    #Output: Returns the new edge length after the conformal varations have been applied
    #Author: MELT, 10/28/14
    #Change Log:
    #def conformalize(self,vertexAVar,vertexBVar,length):  #No longer needed, as conformal variations are always zero
    #    return math.exp(.5*(vertexAVar+vertexBVar))*length
    #

    #input: table of edge objects
    #output: none, debugging function
    #author: MELT, 10/7/2014
    #change log:
    def showEdgeTable(self,tableOfEdges):  #Prints all information stored in edge table, used for depbugging
        count = 0
        for i in range(len(tableOfEdges)):
            for j in range(len(tableOfEdges[i])):
                if tableOfEdges[i][j] != 0:
                    print(tableOfEdges[i][j].vertex1,tableOfEdges[i][j].vertex2,tableOfEdges[i][j].tetrahedraEdgeIsIn, tableOfEdges[i][j].edgelength)
                    count = count+1
        print(count)
    #

    #input: list of tetrahedra names
    #output: number of vertices, builds tetrahedra objects
    #author: TM, 10/3/2014
    #change log: MELT 10/16/14 added calculation to find number of verticies
    def createTetrahedraList(self,primalList):  #creates list of tetrahedra in metric
        numberOfVertices = 0
        for i in range(len(primalList)):
            vertex1 = primalList[i][0]
            vertex2 = primalList[i][1]
            vertex3 = primalList[i][2]
            vertex4 = primalList[i][3]
            self.tetrahedralist.append(Tetrahedron(vertex1,vertex2,vertex3,vertex4))
            self.tetrahedralist[i].edgesintetrahedron.append([vertex1,vertex2])
            self.tetrahedralist[i].edgesintetrahedron.append([vertex1,vertex3])
            self.tetrahedralist[i].edgesintetrahedron.append([vertex1,vertex4])
            self.tetrahedralist[i].edgesintetrahedron.append([vertex2,vertex3])
            self.tetrahedralist[i].edgesintetrahedron.append([vertex2,vertex4])
            self.tetrahedralist[i].edgesintetrahedron.append([vertex3,vertex4])
            largestNumber = self.tetrahedralist[i].vertex4
            if numberOfVertices < largestNumber:
                numberOfVertices = largestNumber
        return numberOfVertices
    #10/16/2014 MELT added vertex counter

    #inpout: list of tetrahedra and a table of edges
    #output: true if there exists an illegal tetrahedron otherwise false, prints Names of tetrahedra and edge information that cause illegal tetrahedra to form
    #Author: M, 10/16/14
    #change log:
    def advancedCheckLegalTetrahedra(self,listOfTetrahedra,tableOfEdges):
        #legal = True
        everIllegal = False
        for i in range(len(listOfTetrahedra)):
            legal = listOfTetrahedra[i].checkLegalTetrahedron(tableOfEdges)
            if legal == False:
                #print('Illegal Tetrahedron: (',listOfTetrahedra[i].vertex1,',',listOfTetrahedra[i].vertex2,',',listOfTetrahedra[i].vertex3,',',listOfTetrahedra[i].vertex4,')')
                everIllegal = True
                break  #delete statement if uncommenting print command
        # if everIllegal == True:
        #     print(":(")
        return everIllegal
    #
    #input: table of edges
    #output: total length of all edges in table
    #Author: ME
    #Change log: 6/16/15, M added first comments to function
    def findTotalEdgeLength(self,tableOfEdges):
        L = 0
        for i in range(1,self.vertexNumber+1):
            for j in range(len(tableOfEdges[0])):  #looks at every entry in edge table
                if tableOfEdges[j][i] != 0:  #0 indicates no edge present, so will only add edge length if edge is present
                    L = L+tableOfEdges[j][i].edgelength
        return L

    #input: target Vertex, table of edges
    #output: curvature of that vertex
    #author: MELT
    #Change log: 6/16/15 added first comments to code
    def calculateVertexCurvature(self,vertex,tableOfEdges):
        vertexCurvature = 0
        for i in range(len(tableOfEdges[0])):  #looks at every edge associated with the given vertex
            if tableOfEdges[vertex][i] != 0:  #checks both possible forms of edge name (aka (i,j) and (j,i))
                vertexCurvature = vertexCurvature+tableOfEdges[vertex][i].edgecurvature  #adds edge curvature as required by vertex curvature definition
            elif tableOfEdges[i][vertex] != 0:
                vertexCurvature = vertexCurvature+tableOfEdges[i][vertex].edgecurvature
        vertexCurvature = vertexCurvature/2  #part of definition of vertex curvature
        return vertexCurvature


    #input: table of edges, error willing to accept to declare LCSC quantity
    #output: true if metric is LCSC, false if metric is not LCSC
    #author: MELT
    #Change Log: M, 6/16/15 added comments to code
    def checkLCSC(self,tableOfEdges,error=.0001):
        LCSC = True
        listOfLengths = []  #edge X will appear in position X and contain a list of all the lengths of edges connected to vertex X
        for i in range(self.vertexNumber+1):  #creates a spot for each vertex
            listOfLengths.append([0])
        for i in range(len(self.edgeList)):  #for each edge
            edge1 = self.edgeList[i][0]
            edge2 = self.edgeList[i][1]
            listOfLengths[edge1].append(self.edgetable[edge1][edge2].edgelength)  #adds edge length to the two verticies are connected
            listOfLengths[edge2].append(self.edgetable[edge1][edge2].edgelength)
       # for i in range(self.vertexNumber):
       #     print(i,self.calculateVertexCurvature(i,tableOfEdges),self.LEHR * L)
        for i in range(1,self.vertexNumber+1):
            if math.fabs(self.calculateVertexCurvature(i,tableOfEdges)-(self.LEHR * sum(listOfLengths[i])/2)) > (error):  #checks if LCSC is met
                #print(math.fabs(self.calculateVertexCurvature(i,tableOfEdges)-(self.LEHR * L)))
                LCSC = False
                break  #stops loop if LCSC is discovered not to exist
        return LCSC


    #Input: Error tolerance allowed to be considered L-Einstein
    #Output: True if L-Einstein requirments met, false otherwise
    def checkLEinstein(self,error=.0001):
        LEinstein = True
        for i in range(len(self.edgeList)):  # for each edge
            edge1 = self.edgeList[i][0]
            edge2 = self.edgeList[i][1]
            curvature=self.edgetable[edge1][edge2].edgecurvature  #find edge curvature
            LEHRl=self.edgetable[edge1][edge2].edgelength  #find edge length
            LEHRl=LEHRl*self.LEHR  #multiply edge length by LEHR
            if math.fabs(curvature-LEHRl)> error:  #check if L-Einstein met
                LEinstein= False
                break  #stops if L-Einstein found not to exist
        return LEinstein

    #inpout: table of edges and number of verticies
    #output: Creates a list such that the n-1 item of the list corresponds with the sum of all the edges incident on edge n
    #Author: EM, 02/04/15
    #change log: 6/16/15 M, added comments

    def getEdgeSums(self,tableOfEdges,numberOfVertices):
        edgeSums = []
        for i in range(1,numberOfVertices+1):  #for each vertex
            sums = 0
            for j in range(1,numberOfVertices+1):  #look at each edge associated with vertex and add its length up
                if tableOfEdges[i][j] != 0:
                    sums = sums + tableOfEdges[i][j].edgelength
                elif tableOfEdges[j][i] != 0:
                    sums = sums + tableOfEdges[j][i].edgelength
            edgeSums.append(sums)
        return edgeSums

    def getEdgeStarOverEdgeSums(self,tableOfEdges,numberOfVertices):
        edgeSumsStar = []
        for i in range(1,numberOfVertices+1):
            sums = 0
            for j in range(1,numberOfVertices+1):
                if tableOfEdges[i][j] != 0:
                    sums = sums + tableOfEdges[i][j].edgelengthStar/tableOfEdges[i][j].edgelength
                elif tableOfEdges[j][i] != 0:
                    sums = sums + tableOfEdges[j][i].edgelengthStar/tableOfEdges[j][i].edgelength
            edgeSumsStar.append(sums)
        return edgeSumsStar

    #input: none
    #output: none, super-mega-function that does EVERYTHING! prints LEHR
    #author: METAL, 10/8/2014
    #change log: Michael 10/16/4
    def __init__(self, backgroundFile = 'backgroundMetric.txt', manifoldFile = 'manifoldExample.txt'):
        self.edgetable = []#A list of lists of edges such that edge x,y is stored in edgetable[x][y]
        self.tetrahedralist = []#A list of tetrahedron objects
        self.edgeList = []#A list of edge names where edgeList[i]=[[a],[b]]
        self.vertexNumber = 0
        self.LEHR = 10000
        self.LCSC = False
        self.LEinstein = False
        self.vertexCurvatureList = []
        self.sumOfEdgesAtVertex = []
        readFile = open(manifoldFile)
        data = readFile.read()        #Prepares file for read in
        data = data.split("facets :=") #Look up strip to remove white space
        data[1] = data[1].strip('[];')
        data[1] = data[1].split('],[')
        tetrahedron = []
        for i in range(0, len(data[1])):   #List comprehensions
            tetrahedron.append(data[1][i])
        for i in range(0, len(tetrahedron)):
            tetrahedron[i] = tetrahedron[i].split(',')
        readFile.close()
        tetrahedron = [[int(i) for i in tetrahedron[j]] for j in range(len(tetrahedron))] #turns tetrahedron from str to int
        #start of calculations, end of file reading (in this function)
        self.vertexNumber = self.createTetrahedraList(tetrahedron)  #builds list of tetrahedra and sets number of verticies in metric
        self.createEdgeTable(self.edgetable,self.vertexNumber)  #builds empty edge table (all zeros)
        self.fillEdgeTable(self.tetrahedralist,self.edgetable,self.edgeList,backgroundFile)  #places edges in edge table (replaces zeros with edges where needed)
        self.sumOfEdgesAtVertex = self.getEdgeSums(self.edgetable,self.vertexNumber)  #finds the sum of edge lengths at each vertex
        illegalTetrahedrons = self.advancedCheckLegalTetrahedra(self.tetrahedralist,self.edgetable)  #checks for legal tetradra
        if illegalTetrahedrons == True:  #if illegal found, further calculations cannot be done and program stops
            self.good = False
        else:
         #   for i in range(len(self.tetrahedralist)): #333 hessian
         #       self.tetrahedralist[i].initalizeTriangles(self.edgetable)  #finds information about triangles within tetrahedra (333 used in hessian calculations only)
            self.totalLength = self.findTotalEdgeLength(self.edgetable)  #gets total edge length of all edges in metric
            for i in range(len(self.tetrahedralist)):  #for each tetrahedra
                self.tetrahedralist[i].calculateDihedralAngles(self.edgetable)  #finds dihedral angles
                #self.tetrahedralist[i].getTetCenDis()  #333 for hessian
            for i in range(len(self.edgeList)):  #for each edge
                edge1 = self.edgeList[i][0]
                edge2 = self.edgeList[i][1]
                self.edgetable[edge1][edge2].calculateEdgeCurvature(self.tetrahedralist)  #finds edge curvature
                #self.edgetable[edge1][edge2].calculateEdgeLengthStar(self.tetrahedralist)  #333 used in hessian
           # self.edgeStarOverEdgeTotal = self.getEdgeStarOverEdgeSums(self.edgetable,self.vertexNumber)  #333 used in hessian
            self.LEHR = self.findLEHR(self.edgeList,self.edgetable)  #finds LEHR
            #self.showEdgeTable(self.edgetable)
            self.LCSC = self.checkLCSC(self.edgetable)  #checks if LCSC and L-Einstein met
            self.LEinstein = self.checkLEinstein()
            self.good = True
    # 10/16/14 Michael; replaced for loop to check legal tetrahedron to advancedcheckLegalTetrahedron command

class metric:
    def calLEHR(self,convar):
        if any(convar) != 0:
            #print("Not all Zero")
            self.vertexCurvatureList = []
            edgeList = self.background.edgeList
            edgeTable = self.setEdgeLengths(convar)
            self.edgeTable = edgeTable
            tetrahedraList = self.background.tetrahedralist
            self.good = not self.advancedCheckLegalTetrahedra(tetrahedraList,edgeTable)
            if self.good == True:
                for i in range(len(tetrahedraList)):
                    tetrahedraList[i].calculateDihedralAngles(edgeTable)
                for i in range(len(edgeList)):
                    edgeTable[edgeList[i][0]][edgeList[i][1]].calculateEdgeCurvature(tetrahedraList)
                for i in range(1,self.background.vertexNumber+1):
                    self.vertexCurvatureList.append(self.calculateVertexCurvature(i,edgeTable))
                result = self.findLEHR(edgeList,edgeTable)
                self.LEHR = result
                self.isLCSC = self.checkLCSC(edgeTable,edgeList)
                self.isLEinstein = self.checkLEinstein(edgeTable,edgeList)
            else:
                result = 1000
        else:
            #print("All Zero")
            self.LEHR = self.background.LEHR
            result = self.LEHR
            self.isLCSC = self.background.LCSC
            self.isLEinstein = self.background.LEinstein
            self.good = self.background.good

        return result



    def checkLCSC(self,tableOfEdges,edgeList,error=.0001):
            LCSC = True
            listOfLengths = []
            for i in range(self.background.vertexNumber+1):
                listOfLengths.append([0])
            for i in range(len(self.background.edgeList)):
                edgeSpot1 = edgeList[i][0]
                edgeSpot2 = edgeList[i][1]
                lengths = tableOfEdges[edgeSpot1][edgeSpot2].edgelength
                listOfLengths[edgeSpot1].append(lengths)
                listOfLengths[edgeSpot2].append(lengths)
           # for i in range(self.vertexNumber):
           #     print(i,self.calculateVertexCurvature(i,tableOfEdges),self.LEHR * L)
            self.LCSCError = 0
            for i in range(1,self.background.vertexNumber+1):
                if self.LCSCError < math.fabs(self.calculateVertexCurvature(i,tableOfEdges)-self.LEHR * sum(listOfLengths[i])/2):
                        self.LCSCError = self.calculateVertexCurvature(i,tableOfEdges)-(self.LEHR * sum(listOfLengths[i])/2)
                if math.fabs(self.calculateVertexCurvature(i,tableOfEdges)-(self.LEHR * sum(listOfLengths[i])/2)) > (error):
                    #print(math.fabs(self.calculateVertexCurvature(i,tableOfEdges)-(self.LEHR * L)))
                    LCSC = False
                   # break
            return LCSC

    def checkLEinstein(self,tableOfEdges,edgeList,error=.0001):
        LEinstein = True
        for i in range(len(self.background.edgeList)):
                edgeSpot1 = edgeList[i][0]
                edgeSpot2 = edgeList[i][1]
        for i in range(len(edgeList)):
            curvature=tableOfEdges[edgeSpot1][edgeSpot2].edgecurvature
            LEHRl=tableOfEdges[edgeSpot1][edgeSpot2].edgelength
            LEHRl=LEHRl*self.LEHR
            if math.fabs(curvature-LEHRl)> error:
                LEinstein= False
        return LEinstein


    def optimizeLEHR(self,convar):
        res = minimize(self.calLEHR, convar ,method = 'nelder-mead',options={'disp':False})
        #res = minimize(self.calLEHR, convar ,method = 'Newton-CG',jac=self.grad,hess = self.hess,options={'disp':True})
        #res = minimize(self.simplerFunction, self.x0 ,method = 'Newton-CG',jac = simplerFunctionDer,hessp = simplerFunctionHes,options={'disp':True})
        self.minima = res.x
        return res

    #not optimized
    def grad(self,convar):
        Grad = []
        self.calLEHR(convar)
        for i in range(self.background.vertexNumber):
                #creating the gradient
                temp = self.vertexCurvatureList[i]
                temp = temp-self.LEHR*.5*self.sumOfEdgesAtVertexList[i]
                temp = temp/self.totalEdgeLength
                Grad.append(temp)
        return numpy.array(Grad)

    #not optimized, or any hessian related functions
    def hess(self, convar):
        self.vertexCurvatureList = []
        edgeList = copy.deepcopy(self.background.edgeList)
        edgeTable = self.setEdgeLengths(convar)
        tetrahedraList = copy.deepcopy(self.background.tetrahedralist)
        self.good = not self.advancedCheckLegalTetrahedra(tetrahedraList,edgeTable)
        if self.good == True:
            for i in range(len(tetrahedraList)):
                tetrahedraList[i].initalizeTriangles(edgeTable)
                tetrahedraList[i].calculateDihedralAngles(edgeTable)
                tetrahedraList[i].getTetCenDis()
            for i in range(len(edgeList)):
                edgeTable[edgeList[i][0]][edgeList[i][1]].calculateEdgeCurvature(tetrahedraList)
                edgeTable[edgeList[i][0]][edgeList[i][1]].calculateEdgeLengthStar(tetrahedraList)
            for i in range(1,self.background.vertexNumber+1):
                self.vertexCurvatureList.append(self.calculateVertexCurvature(i,edgeTable))
            edgeStarOverEdgeTotal = self.getEdgeStarOverEdgeSums(edgeTable,self.background.vertexNumber)
            result = self.findLEHR(edgeList,edgeTable)
            self.LEHR = result
            self.isLCSC = self.checkLCSC(edgeTable,edgeList)
        else:
            result = 1000
        hessian = []
        for i in range(self.background.vertexNumber):
        #creating the hessian
            hessian.append([0]*self.background.vertexNumber)
            for j in range(self.background.vertexNumber):
                if edgeTable[i+1][j+1] == 0 and edgeTable[j+1][i+1] == 0 and i!=j:
                    if edgeTable[i+1][j+1] == 0:
                        hessian[i][j] = 0
                elif i == j:
                    hessian[i][j] = self.hessianSame(edgeStarOverEdgeTotal,i+1)
                else:
                    hessian[i][j] = self.hessianDifferent(edgeTable,i+1,j+1)
        self.makePosDef(numpy.asarray(hessian),.1)
        return hessian

    def hessianDifferent(self,edgetable,vertex1,vertex2):
        if vertex1 > vertex2:
            lijStarOverlij = edgetable[vertex2][vertex1].edgelengthStar/edgetable[vertex2][vertex1].edgelength
            edgeCurvature = edgetable[vertex2][vertex1].edgecurvature
            edgeLength = edgetable[vertex2][vertex1].edgelength
        else:
            lijStarOverlij = edgetable[vertex1][vertex2].edgelengthStar/edgetable[vertex1][vertex2].edgelength
            edgeCurvature = edgetable[vertex1][vertex2].edgecurvature
            edgeLength = edgetable[vertex1][vertex2].edgelength
        partial = -2*lijStarOverlij/self.totalEdgeLength
        partial = partial + edgeCurvature/(4*self.totalEdgeLength)
        partial = partial - self.LEHR*edgeLength/(4*self.totalEdgeLength)
        partial = partial - self.vertexCurvatureList[vertex2-1]*self.sumOfEdgesAtVertexList[vertex1-1]/(2*self.totalEdgeLength**2)
        partial = partial + self.sumOfEdgesAtVertexList[vertex1-1]*self.LEHR*self.sumOfEdgesAtVertexList[vertex2-1]/(2*self.totalEdgeLength**2)
        partial = partial - self.sumOfEdgesAtVertexList[vertex2-1]*self.vertexCurvatureList[vertex1-1]/(2*self.totalEdgeLength**2)
        return partial

    def hessianSame(self,edgeStarOverEdgeTotal,vertex):
        partial = 2*edgeStarOverEdgeTotal[vertex-1]/self.totalEdgeLength
        partial = partial + self.vertexCurvatureList[vertex-1]/(2*self.totalEdgeLength)
        partial = partial - self.LEHR*self.sumOfEdgesAtVertexList[vertex-1]**2/(4*self.totalEdgeLength)
        partial = partial - self.sumOfEdgesAtVertexList[vertex-1]*self.vertexCurvatureList[vertex-1]/(self.totalEdgeLength**2)
        partial = partial + self.LEHR*self.sumOfEdgesAtVertexList[vertex-1]**2/(2*self.totalEdgeLength**2)
        return partial

    def makePosDef(self,target,add = 1):
        eigns = numpy.linalg.eigvals((target))
        smallEig = min(eigns)
        if(smallEig <= 0):
            target = target + (math.fabs(smallEig)+add)*numpy.identity(target.size**.5)
        eigns2=numpy.linalg.eigvals(target)
        return target

    #true if there exists an illegal tetrahedron
    def advancedCheckLegalTetrahedra(self,listOfTetrahedra,tableOfEdges):
            #legal = True
            everIllegal = False
            for i in range(len(listOfTetrahedra)):
                legal = listOfTetrahedra[i].checkLegalTetrahedron(tableOfEdges)
                if legal == False:
                    #print('Illegal Tetrahedron: (',listOfTetrahedra[i].vertex1,',',listOfTetrahedra[i].vertex2,',',listOfTetrahedra[i].vertex3,',',listOfTetrahedra[i].vertex4,')')
                    everIllegal = True
                    break  #comment this out if print commands are used
            return everIllegal

    def findLEHR(self,listOfEdges,tableOfEdges):
            listOfLengths = []
            listOfCurvatures = []
            for i in range(len(listOfEdges)):
                edgeSpot1 = listOfEdges[i][0]
                edgeSpot2 = listOfEdges[i][1]
                listOfLengths.append(tableOfEdges[edgeSpot1][edgeSpot2].edgelength)
                listOfCurvatures.append(tableOfEdges[edgeSpot1][edgeSpot2].edgecurvature)
            return sum(listOfCurvatures)/sum(listOfLengths)

    def calculateVertexCurvature(self,vertex,tableOfEdges):
        vertexCurvature = 0
        for i in range(len(tableOfEdges[0])):
            if tableOfEdges[vertex][i] != 0:
                vertexCurvature = vertexCurvature+tableOfEdges[vertex][i].edgecurvature
            elif tableOfEdges[i][vertex] != 0:
                vertexCurvature = vertexCurvature+tableOfEdges[i][vertex].edgecurvature
        vertexCurvature = vertexCurvature/2
        return vertexCurvature

    def setEdgeLengths(self,convar):
        edgeLengthTable = copy.deepcopy(self.background.edgetable)
        for i in range(len(edgeLengthTable)):
            for j in range(len(edgeLengthTable[i])):
                if(edgeLengthTable[i][j]!=0):
                    length = edgeLengthTable[i][j].edgelength
                    lengthOverTwo = length/2
                    edgeLengthTable[i][j].edgelength = math.exp(.5*(convar[i-1]+convar[j-1]))*length
                    self.totalEdgeLength = self.totalEdgeLength + length/2
                    self.sumOfEdgesAtVertexList[i-1] = self.sumOfEdgesAtVertexList[i-1]+lengthOverTwo
                    self.sumOfEdgesAtVertexList[j-1] = self.sumOfEdgesAtVertexList[j-1]+lengthOverTwo
        return edgeLengthTable

    def getEdgeStarOverEdgeSums(self,tableOfEdges,numberOfVertices):
            edgeSumsStar = []
            for i in range(1,numberOfVertices+1):
                sums = 0
                for j in range(1,numberOfVertices+1):
                    if tableOfEdges[i][j] != 0:
                        sums = sums + tableOfEdges[i][j].edgelengthStar/tableOfEdges[i][j].edgelength
                        #print(str(i)+","+str(j)+": "+str(tableOfEdges[i][j].edgelengthStar))
                    elif tableOfEdges[j][i] != 0:
                        sums = sums + tableOfEdges[j][i].edgelengthStar/tableOfEdges[j][i].edgelength
                        #print(str(j)+","+str(i)+": "+str(tableOfEdges[j][i].edgelengthStar))
                edgeSumsStar.append(sums)
            return edgeSumsStar

    def __init__(self,backgroundFile,triagulation,convar):
        self.background = backgroundMetricClass(backgroundFile,triagulation)
        self.totalEdgeLength = 0
        self.sumOfEdgesAtVertexList = [0]*self.background.vertexNumber
        self.vertexCurvatureList = [0]
        self.isLCSC = False
        self.isLEinstein = False
        self.good = True
        self.LEHR = 1000
        self.edgeTable = []
        self.LCSCError = 999
        self.minima = []


def ToReducedRowEchelonForm( M):
    if not M: return
    lead = 0
    rowCount = len(M)
    columnCount = len(M[0])
    for r in range(rowCount):
        if lead >= columnCount:
            return
        i = r
        while M[i][lead] == 0:
            i += 1
            if i == rowCount:
                i = r
                lead += 1
                if columnCount == lead:
                    return
        M[i],M[r] = M[r],M[i]
        lv = M[r][lead]
        M[r] = [ mrx / float(lv) for mrx in M[r]]
        for i in range(rowCount):
            if i != r:
                lv = M[i][lead]
                M[i] = [ iv - lv*rv for rv,iv in zip(M[r],M[i])]
        lead += 1


def modifyBackground(cList,base,filename):
    backgroundMetric = open(filename,"r")
    storage = backgroundMetric.readlines()
    nameList = []
    lengthList = []
    i = 0
    while(i < len(storage)-1):
        #breaks background metric into a names and lengths, such that name at position i corresponds to length and position i
        nameList.append(storage[i])
        lengthList.append(float(storage[i+1]))
        i = i+2
    # makes files usable
    for i in range(len(nameList)):
        nameList[i] = nameList[i].split(",")
        nameList[i][0] = int(nameList[i][0])
        nameList[i][1] = int(nameList[i][1])
        #nameList[i].append(0)  #Possible place for error, does nameList[i]  correspond to basis edge[i]
    base = numpy.array(base)
    cList = numpy.array(cList)
    lengthList = numpy.dot(base.transpose(),cList)
    for i in range(len(lengthList)):
        lengthList[i] = math.exp(lengthList[i])**.5
    #lengthList = [(math.exp(lengthList[i]))**.5 for i in lengthList]
#    lengthList = [(math.exp(c1-c4+c5))**.5,(math.exp(c1))**.5,(math.exp(-c2))**.5,(math.exp(-c1+c4-c5))**.5,(math.exp(-c2+c3+c4))**.5,(math.exp(-c5))**.5,(math.exp(-c3))**.5,(math.exp(-c1+c2-c3))**.5,(math.exp(-c4))**.5,(math.exp(c1+c3+c5))**.5]
    backgroundMetric = open(filename,"w")
    for i in range(len(nameList)):
        backgroundMetric.write(str(nameList[i][0])+","+str(nameList[i][1])+"\n"+str(lengthList[i])+"\n")
       # print(str(nameList[i][0])+","+str(nameList[i][1])+"\n"+str(lengthList[i])+"\n")
    backgroundMetric.close()

def generateConvar(size,numberSame):
    numberSame = numberSame%(size+1)
    conVar = []
    temp = []
    for i in range(size-numberSame):
        conVar.append(random.randint(-1,3))
    temp = random.randint(-1,3)
    for i in range(numberSame):
        conVar.append(temp)
    conVar[0] = random.random()
    random.shuffle(conVar)
    return conVar

def advancedGenerateConvar(size,background,triagulation,givenConVar = [], steps = 17, startingSize = 1,outerLoop = 1):
    findImprove = False
    if len(givenConVar) == 0:
        conVarZero = [0] * size
        finalConVar = [0] * size
    else:
        conVarZero = copy.deepcopy(givenConVar)
        finalConVar = copy.deepcopy(givenConVar)
    storage = 0
    maxTemp = outerLoop
    for k in range(outerLoop):
        for i in range(size):
            stepSize = startingSize
            storage = conVarZero[i]
            for j in range(steps):
                base = metric(background,triagulation,conVarZero)
                base.calLEHR(conVarZero)
                conVarZero[i] = conVarZero[i]+stepSize
                temp1 = metric(background,triagulation,conVarZero)
                temp1.calLEHR(conVarZero)
                conVarZero[i] = conVarZero[i] - 2 * stepSize
                temp2 = metric(background,triagulation,conVarZero)
                temp2.calLEHR(conVarZero)
                if min(temp1.LEHR,temp2.LEHR) > base.LEHR:
                    conVarZero[i] = conVarZero[i] + stepSize
                elif temp1.LEHR < temp2.LEHR:
                    conVarZero[i] = conVarZero[i] + 2 * stepSize
                    findImprove = True
                else:
                    findImprove = True
                if findImprove == False:
                    stepSize = stepSize/2
                    findImprove = False
                print(base.LEHR)
            finalConVar[i] = conVarZero[i]
            #conVarZero[i] = storage
        base = metric(background,triagulation,finalConVar)
        base.calLEHR(finalConVar)
        conVarZero = copy.deepcopy(finalConVar)
    print(base.LEHR)
    print("\n")
    return  finalConVar


def pentachoronWalk(basis,numberVerts,backgroundfile,triangulation,restarts = 100,numberBackgrounds = 10):
    results = []
    failures = []
    times = []
    lengths = []
    volumes = []
    failedlengths =[]
    working = True
    temp = []
    for i in range(len(basis)):
        temp.append('c'+str(i+1))
    for i in range(numberVerts):
        temp.append('f'+str(i+1))
    temp.append("LEHR")
    temp.append("LCSC")
    temp.append("L-Einstein")
    temp.append("numberRestarts")
    temp.append("LCSC Error")
    results.append(temp)
    temp = []
    for i in range(len(basis)):
        temp.append('c'+str(i+1))
    for i in range(numberVerts):
        temp.append('f'+str(i+1))
    temp.append("LEHR")
    temp.append("LCSC")
    temp.append("Is LEinstein")
    temp.append("LCSC Error")
    failures.append(temp)
    times.append(["Full Search Time","Generate Moduli Time","Optomize Time","Starting State/Conformal Varition Finder Time","Restarts"])
    lengths.append(["Edges"])
    startAll = time.time()
    for k in range(numberBackgrounds):
        walkThroughBack(basis,numberVerts,backgroundfile,triangulation,numberBackgrounds,results,times,lengths,failures,volumes,failedlengths,k,restarts)
    endAll = time.time()
    print(endAll-startAll)
    return [results,failures,times,lengths,volumes,failedlengths]

def legalBackground(cList,base,background,triagulation):
    modifyBackground(cList,base,background)
    cTest = [0]*len(cList)
    test = metric(background,triagulation,cTest)
    LEHR = test.calLEHR(cTest)
    if LEHR >= 1000:
       return [False,0]
    else:
        return [True,test]
def walkThroughBack(basis,numberVerts,backgroundfile,triangulation,numberBackgrounds,results,times,lengths,failures,volumes,failedlengths,k,restarts):
    startBackground = time.time()
    happyBackground = [False,0]
    startModuli = time.time()
    working = True
    happyConVar = False
    startConar = time.time()
    cList = []
    for i in range(len(basis)):
        cList.append(random.random()*20-1)
    squares = [i**2 for i in cList]
    while sum(squares) > 10**2 or happyBackground[0] == False:
        # if sum(squares) > 10**2:
        #     print("Finding C-List")
        # else:
        #     print("NOt")
        cList = []
        for i in range(len(basis)):
            cList.append(random.random()-.5)
        squares = [i**2 for i in cList]
        happyBackground = legalBackground(cList, basis,backgroundfile,triangulation)
    endModuli = time.time()
    conVar = [0]*numberVerts
    print("background number "+str(k+1) +" out of " + str(numberBackgrounds))
    startOpt = time.time()
    #add optimize
    test = happyBackground[1]
    temp = test.optimizeLEHR(conVar)
    conVar = temp.x
    test.calLEHR(conVar)
    endOpt = time.time()
    working = test.isLCSC
    if working == True:
        listOfLengths = []
        conVar = list(conVar)
        results.append([cList,conVar,test.LEHR,test.isLCSC,test.isLEinstein,-1,test.LCSCError])
        endBackground = time.time()
        times.append([endBackground-startBackground,endModuli-startModuli,endOpt-startOpt,0,0])
        for i in range(len(test.edgeTable)):
            for j in range(len(test.edgeTable)):
                if test.edgeTable[i][j] != 0:
                    listOfLengths.append(test.edgeTable[i][j].edgelength)
        lengths.append(listOfLengths)
    else:
        for j in range(restarts):
            happyConVar = False
            same = j
            happyConVar = test.good
            while(happyConVar == False):
                conVar = generateConvar(numberVerts,same)
                test = metric(backgroundfile,triangulation,conVar)
                test.calLEHR(conVar)
                happyConVar = test.good
            conVar = advancedGenerateConvar(numberVerts,backgroundfile,triangulation,conVar)
            print(conVar)
            test = metric(backgroundfile,triangulation,conVar)
            test.calLEHR(conVar)
            endConVar = time.time()
            startOpt = time.time()
            temp = test.optimizeLEHR(conVar)
#            test.calLEHR(conVar)  #Important?
            endOpt = time.time()
            conVar = temp.x
            test.calLEHR(conVar)
            working = test.isLCSC
            if working == True:
                listOfLengths = []
                conVar = list(conVar)
                results.append([cList,conVar,test.LEHR,test.isLCSC,test.isLEinstein,j,test.LCSCError])
                endBackground = time.time()
                times.append([endBackground-startBackground,endModuli-startModuli,endOpt-startOpt,endConVar-startConar,j])
                for i in range(len(test.edgeTable)):
                    for j in range(len(test.edgeTable)):
                        if test.edgeTable[i][j] != 0:
                            listOfLengths.append(test.edgeTable[i][j].edgelength)
                lengths.append(listOfLengths)
                # listOfVolumes = [] #uncomment block to get the volumes of metrics that find L-CSC
                # for i in range(len(test.background.tetrahedralist)):
                #     vol = test.background.tetrahedralist[i].volumeOfTetrahedron(test.edgeTable)
                #     listOfVolumes.append(vol)
                # volumes.append(listOfVolumes)
                break
    if working == False:
        conVar = list(conVar)
        failures.append([cList,conVar,test.LEHR,test.isLCSC,test.isLEinstein,test.LCSCError])
        endBackground = time.time()
        times.append([endBackground-startBackground,endModuli-startModuli,endOpt-startOpt])
        listOfVolumes = []
        for i in range(len(test.background.tetrahedralist)):
            vol = test.background.tetrahedralist[i].volumeOfTetrahedron(test.edgeTable)
            listOfVolumes.append(vol)
        volumes.append(listOfVolumes)
        listOfFailedLenghts = []
        for i in range(len(test.edgeTable)):
            for j in range(len(test.edgeTable)):
                if test.edgeTable[i][j] != 0:
                    listOfFailedLenghts.append(test.edgeTable[i][j].edgelength)
        failedlengths.append(listOfFailedLenghts)
        #print(temp)
def main():
    start= time.time()
    storage = str(0)+".txt"
    LEHRList = []
    numberVertices=5
    numberOfBackgrounds=10
    numberRestarts = 1
    #seed=4741252
    seed = 32190
    #seed=263594
    #seed=56932684
    #seed=71293
    random.seed(seed)
    #seed=9865721
    triangulation='manifoldExample4.txt'
    backgroundfile='backgroundMetric.txt'
    basis = BasisBuilder.main(triangulation)
    #basis = numpy.array([[1,1,0,-1,0,0,0,-1,0,1],[0,0,-1,0,-1,0,0,1,0,0],[0,0,0,0,1,0,-1,-1,0,1],[-1,0,0,1,1,0,0,0,-1,0],[1,0,0,-1,0,-1,0,0,0,1]])
    faceInfo = " "
    print("Hello World!\n")
    results = pentachoronWalk(basis,numberVertices,backgroundfile,triangulation,numberRestarts,numberOfBackgrounds)
    for i in range(len(results[0])):
        print(results[0][i])
    for i in range(len(results[1])):
        print(results[1][i])
    foundResultsFile = open("foundResults.txt","w")
    for i in range(1,len(results[0])):
        for j in range(len(results[0][i])):
            foundResultsFile.write(str(results[0][i][j])+" ")
        foundResultsFile.write("\n")
    foundResultsFile.close()
    notfoundResultsFile = open("notFoundResults.txt","w")
    for i in range(1,len(results[1])):
        for j in range(len(results[1][i])):
            notfoundResultsFile.write(str(results[1][i][j])+" ")
        notfoundResultsFile.write("\n")
    notfoundResultsFile.close()
    lengthsFile = open("foundLengths.txt","w")
    for i in range(1,len(results[3])):
        for j in range(len(results[3][i])):
            lengthsFile.write(str(results[3][i][j])+" ")
        lengthsFile.write("\n")
    lengthsFile.close()
    volumeFile = open("NotFoundVolumes.txt","w")
    for i in range(0,len(results[4])):
        for j in range(len(results[4][i])):
            volumeFile.write(str(results[4][i][j])+" ")
        volumeFile.write("\n")
    volumeFile.close()
    notFoundLengthsFile = open("notFoundLengths.txt","w")
    for i in range(0,len(results[5])):
        for j in range(len(results[5][i])):
            notFoundLengthsFile.write(str(results[5][i][j])+" ")
        notFoundLengthsFile.write("\n")
    notFoundLengthsFile.close()
    timeFile = open("timeInformation.txt","w")
    for i in range(1,len(results[2])):
        for j in range(len(results[2][i])):
            timeFile.write(str(results[2][i][j])+" ")
        timeFile.write("\n")
    end=time.time()
    print(end-start)
    timeFile.write("\n\n"+str(end-start))
    timeFile.close()
    notfoundResultsFile.close()

main()
#cProfile.run('main()')