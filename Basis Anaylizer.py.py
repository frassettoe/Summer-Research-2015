__author__ = 'Owner'


#tp://www.ripon.edu/wp-content/uploads/2012/10/doubletetra.pdf
#10/16/14 Michael: Added advancedCheckLegalTetrahedra function for checking legal tetrahedron and used it main, made minor change checkLegalTetrahedron
import math
import random
import copy
import numpy

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
        self.tetrahedraEdgeIsIn = [] #This is a number that corresponds to the location of the tetrahedron in a list of Tetrahedron
        self.edgelengthStar = 0
    #10/3/2014: added tetrahedraEdgeIsIn, MELT
    #10/27/14: Vertex class now takes a third argument, length, with default value of one
    #this is just a test ignore this comment
    #Now testing the different branches


    #input: list of Tetrahedron objects
    #output: none, sets edgecurvature to the curvature of the edge
    #author: MELT, 10/14/2014
    #change log: none
    def calculateEdgeCurvature(self,listOfTetrahedra):
        diList = []
        for i in range (len(self.tetrahedraEdgeIsIn)):
            tetLocation = self.tetrahedraEdgeIsIn[i]
            # Finds the dihedral angle of edge in tetrahedron at tetrahedraEdgeIsIn[i]
            #Dihedral angle N at position L corresponds with position of Edge N at position L in edgesintetrahedron
            singleDi = listOfTetrahedra[tetLocation].dihedralanglelist[listOfTetrahedra[tetLocation].edgesintetrahedron.index([self.vertex1,self.vertex2])]
            diList.append(singleDi)
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
            hijk = listOfTetrahedra[tetLocation].faceList[listOfTetrahedra[tetLocation].vertexList.index(lLoc)].hList[listOfTetrahedra[tetLocation].faceList[listOfTetrahedra[tetLocation].vertexList.index(lLoc)].vertexList.index(kLoc)]
            hijkl = listOfTetrahedra[tetLocation].tetCenDisList[listOfTetrahedra[tetLocation].vertexList.index(lLoc)]
            hijl = listOfTetrahedra[tetLocation].faceList[listOfTetrahedra[tetLocation].vertexList.index(kLoc)].hList[listOfTetrahedra[tetLocation].faceList[listOfTetrahedra[tetLocation].vertexList.index(kLoc)].vertexList.index(lLoc)]
            hijlk = listOfTetrahedra[tetLocation].tetCenDisList[listOfTetrahedra[tetLocation].vertexList.index(kLoc)]
            singleStar = hijk*hijkl+hijl*hijlk
            starList.append(singleStar)
        self.edgelengthStar = .5*sum(starList)


class face:
    def getH(self,leg1,leg2,target):
        print('Hi')

    def getAngle(self,c,a,b):
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
        self.angleList = [math.acos(self.getAngle(self.edgeLength[0],self.edgeLength[1],self.edgeLength[2])),math.acos(self.getAngle(self.edgeLength[1],self.edgeLength[0],self.edgeLength[2])),math.acos(self.getAngle(self.edgeLength[2],self.edgeLength[1],self.edgeLength[0]))]
        self.hList = [self.getTriCenDis(self.edgeLength[2],self.edgeLength[0],self.angleList[1]),self.getTriCenDis(self.edgeLength[0],self.edgeLength[1],self.angleList[2]),self.getTriCenDis(self.edgeLength[1],self.edgeLength[2],self.angleList[0])]


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
        self.edgesintetrahedron = []  #stored as lists ex, edge 1,2 stored as [1,2]
        self.dihedralanglelist = []
        self.tetCenDisList = [0]*4
        self.faceList = [0]*4 #faceList[i] is the face oppositve vertex i-1
    #


    def initalizeTriangles(self,tableOfEdges):
        self.faceList[0] = face(tableOfEdges,self.vertex2,self.vertex3,self.vertex4)
        self.faceList[1] = face(tableOfEdges,self.vertex1,self.vertex3,self.vertex4)
        self.faceList[2] = face(tableOfEdges,self.vertex1,self.vertex2,self.vertex4)
        self.faceList[3] = face(tableOfEdges,self.vertex1,self.vertex2,self.vertex3)


    #input: Table of edge objects
    #output: finalDecision, a boolean that tells whether legality conditions have been met
    #author: MELT, 10/13/2014
    #change log: Michael 10/16/14
    def checkLegalTetrahedron(self,tableOfEdges):
        legalTriangle=True
        edge1=tableOfEdges[self.vertex1][self.vertex2].edgelength
        edge2=tableOfEdges[self.vertex1][self.vertex3].edgelength
        edge3=tableOfEdges[self.vertex1][self.vertex4].edgelength
        edge4=tableOfEdges[self.vertex2][self.vertex3].edgelength
        edge5=tableOfEdges[self.vertex2][self.vertex4].edgelength
        edge6=tableOfEdges[self.vertex3][self.vertex4].edgelength
        if edge1+edge2+edge4-2*max(edge1,edge2,edge4) <= 0:
            legalTriangle=False
        if edge1+edge3+edge5-2*max(edge1,edge3,edge5) <= 0:
            legalTriangle=False
        if edge2+edge3+edge6-2*max(edge2,edge3,edge6) <= 0:
            legalTriangle=False
        if edge4+edge5+edge6-2*max(edge4,edge5,edge6) <= 0:
            legalTriangle=False
        legalTetrahedron=True
        det = (-2)*(edge1**4)*(edge6**2)-(2)*(edge1**2)*(edge2**2)*(edge4**2)+(2)*(edge1**2)*(edge2**2)*(edge5**2)+(2)*(edge1**2)*(edge2**2)*(edge6**2)+(2)*(edge1**2)*(edge3**2)*(edge4**2)-(2)*(edge1**2)*(edge3**2)*(edge5**2)+(2)*(edge1**2)*(edge3**2)*(edge6**2)+(2)*(edge1**2)*(edge4**2)*(edge6**2)+(2)*(edge1**2)*(edge5**2)*(edge6**2)-(2)*(edge1**2)*(edge6**4)-(2)*(edge2**4)*(edge5**2)+(2)*(edge2**2)*(edge3**2)*(edge4**2)+(2)*(edge2**2)*(edge3**2)*(edge5**2)-(2)*(edge2**2)*(edge3**2)*(edge6**2)+(2)*(edge2**2)*(edge4**2)*(edge5**2)-(2)*(edge2**2)*(edge5**4)+(2)*(edge2**2)*(edge5**2)*(edge6**2)-(2)*(edge3**4)*(edge4**2)-(2)*(edge3**2)*(edge4**4)+(2)*(edge3**2)*(edge4**2)*(edge5**2)+(2)*(edge3**2)*(edge4**2)*(edge6**2)-(2)*(edge4**2)*(edge5**2)*(edge6**2)
        if det<=0:
            legalTetrahedron=False
        finalDecision=legalTetrahedron and legalTriangle
        return finalDecision
    # 10/16/14; Michael; change return(finalDecision) to return finalDecision

    #input: c,a,b, which are edge lengths
    #output: temp, the cosine of the angle opposite of edge c
    #author: MELT, 10/9/2014
    #change log:
    def getAngle(self,c,a,b):
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
    def calDiAngle(self,Eij,Eik,Eil,Ejk,Ejl,Ekl):
        temp = self.getAngle(Ekl,Eik,Eil)-self.getAngle(Ejk,Eij,Eik)*self.getAngle(Ejl,Eij,Eil)
        temp = temp/math.sin(math.acos(self.getAngle(Ejk,Eij,Eik)))
        temp = temp/math.sin(math.acos(self.getAngle(Ejl,Eij,Eil)))
        return temp
    #

    #input: table of edge objects
    #output: none, finds the dihedral angle at each edge and adds it to the dihedralanglelist
    #author: MELT, 10/9/2014
    #change log:
    def calculateDihedralAngles(self,tableOfEdges):
        edge12=tableOfEdges[self.vertex1][self.vertex2].edgelength
        edge13=tableOfEdges[self.vertex1][self.vertex3].edgelength
        edge14=tableOfEdges[self.vertex1][self.vertex4].edgelength
        edge23=tableOfEdges[self.vertex2][self.vertex3].edgelength
        edge24=tableOfEdges[self.vertex2][self.vertex4].edgelength
        edge34=tableOfEdges[self.vertex3][self.vertex4].edgelength
        self.dihedralanglelist.append(math.acos(self.calDiAngle(edge12,edge13,edge14,edge23,edge24,edge34)))
        self.dihedralanglelist.append(math.acos(self.calDiAngle(edge13,edge12,edge14,edge23,edge34,edge24)))
        self.dihedralanglelist.append(math.acos(self.calDiAngle(edge14,edge13,edge12,edge34,edge24,edge23)))
        self.dihedralanglelist.append(math.acos(self.calDiAngle(edge23,edge12,edge24,edge13,edge34,edge14)))
        self.dihedralanglelist.append(math.acos(self.calDiAngle(edge24,edge12,edge23,edge14,edge34,edge13)))
        self.dihedralanglelist.append(math.acos(self.calDiAngle(edge34,edge13,edge23,edge14,edge24,edge12)))
        # print(self.edgesintetrahedron[0],self.dihedralanglelist[0],self.vertex1,self.vertex2)
        # print(self.edgesintetrahedron[1],self.dihedralanglelist[1],self.vertex1,self.vertex3)
        # print(self.edgesintetrahedron[2],self.dihedralanglelist[2],self.vertex1,self.vertex4)
        # print(self.edgesintetrahedron[3],self.dihedralanglelist[3],self.vertex2,self.vertex3)
        # print(self.edgesintetrahedron[4],self.dihedralanglelist[4],self.vertex2,self.vertex4)
        # print(self.edgesintetrahedron[5],self.dihedralanglelist[5],self.vertex3,self.vertex4)
    #



class metric:
    #input: list of edge names and table of edge objects
    #output: a number that is LEHR
    #author: MELT, 10/9/2014
    #change log:
    def findLEHR(self,listOfEdges,tableOfEdges):
        listOfLengths = []
        listOfCurvatures = []
        for i in range(len(listOfEdges)):
            listOfLengths.append(tableOfEdges[listOfEdges[i][0]][listOfEdges[i][1]].edgelength)
            listOfCurvatures.append(tableOfEdges[listOfEdges[i][0]][listOfEdges[i][1]].edgecurvature)
        return sum(listOfCurvatures)/sum(listOfLengths)
    #

    #input: empty table of edges, number of vertices
    #output: none, makes the edge table of proper size, i.e. 15+1
    #author: MELT, 10/3/14
    #change log: M 10/3/14
    def createEdgeTable(self,edgetable,numberOfVertices = 15):
        numberOfVertices = numberOfVertices+1
        for row in range(numberOfVertices):
            edgetable.append([])
            for column in range(numberOfVertices):
                edgetable[row].append(0)
    #10/3/14 M, edgetable is now an input

    #input: list of tetrahedron object, table of edge objects, list of edge names, file of background metric created from LEHRBackgroundMetricBuilder, word NONE assumes all edge lengths one and no file is real
    #output: none, fills the edge table with edge objects and puts the tetrahedron into the list of tetrahedraEdgeIsIn
    #author: MELT, 10/6/2014
    #change log: 11/18/14 ME
    def fillEdgeTable(self,listOfTetrahedra,tableOfEdges,listOfEdges,fileName = "backgroundMetric.txt",conformalVariationList = "backgroundMetric.txt"):
        #Reads in a background metric file
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
        # makes files usable
        for i in range(len(nameList)):
            nameList[i] = nameList[i].split(",")
            nameList[i][0] = int(nameList[i][0])
            nameList[i][1] = int(nameList[i][1])
        counter = 0
        for i in range(len(listOfTetrahedra)):
            for j in range(len(listOfTetrahedra[i].edgesintetrahedron)):
                    # If edge is not in the the tableOfEdges, add it
                    if tableOfEdges[listOfTetrahedra[i].edgesintetrahedron[j][0]][listOfTetrahedra[i].edgesintetrahedron[j][1]] == 0:
                        # Assigns edge its name
                        tableOfEdges[listOfTetrahedra[i].edgesintetrahedron[j][0]][listOfTetrahedra[i].edgesintetrahedron[j][1]] = Edge(listOfTetrahedra[i].edgesintetrahedron[j][0],listOfTetrahedra[i].edgesintetrahedron[j][1],self.conformalize(conformalVariationList[nameList[counter][0]-1],conformalVariationList[nameList[counter][1]-1],lengthList[counter]))
                        counter = counter+1
                        # Adds tetrahedran to list of tetrahedra edge is in
                        tableOfEdges[listOfTetrahedra[i].edgesintetrahedron[j][0]][listOfTetrahedra[i].edgesintetrahedron[j][1]].tetrahedraEdgeIsIn.append(i)
                        listOfEdges.append([listOfTetrahedra[i].edgesintetrahedron[j][0],listOfTetrahedra[i].edgesintetrahedron[j][1]])
                    # adds edge to tetrahera list if edge already exists
                    else:
                        tableOfEdges[listOfTetrahedra[i].edgesintetrahedron[j][0]][listOfTetrahedra[i].edgesintetrahedron[j][1]].tetrahedraEdgeIsIn.append(i)
    #10/29/14: Michael Added new ability! Function can now be given a backgroundmetric and a list of conformal variations.
    #Both must be given in order for new function attibutes to be used.  Old use without conformal variations and background metric can still be used.
    #11/18/14: Micahel and Erin removed the ability to not take a background metric and conformal variations.

    #Input: The two vertex conformal varations and the length of the edge in question
    #Output: Returns the new edge length after the conformal varations have been applied
    #Author: MELT, 10/28/14
    #Change Log:
    def conformalize(self,vertexAVar,vertexBVar,length):
        return math.exp(.5*(vertexAVar+vertexBVar))*length
    #

    #input: table of edge objects
    #output: none, debugging function
    #author: MELT, 10/7/2014
    #change log:
    def showEdgeTable(self,tableOfEdges):
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
    def createTetrahedraList(self,primalList):
        numberOfVertices = 0
        for i in range(len(primalList)):
            self.tetrahedralist.append(Tetrahedron(primalList[i][0],primalList[i][1],primalList[i][2],primalList[i][3]))
            self.tetrahedralist[i].edgesintetrahedron.append([primalList[i][0],primalList[i][1]])
            self.tetrahedralist[i].edgesintetrahedron.append([primalList[i][0],primalList[i][2]])
            self.tetrahedralist[i].edgesintetrahedron.append([primalList[i][0],primalList[i][3]])
            self.tetrahedralist[i].edgesintetrahedron.append([primalList[i][1],primalList[i][2]])
            self.tetrahedralist[i].edgesintetrahedron.append([primalList[i][1],primalList[i][3]])
            self.tetrahedralist[i].edgesintetrahedron.append([primalList[i][2],primalList[i][3]])
            if numberOfVertices < self.tetrahedralist[i].vertex4:
                numberOfVertices = self.tetrahedralist[i].vertex4
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
        return everIllegal
    #

    def findTotalEdgeLength(self,tableOfEdges):
        L = 0
        for i in range(1,self.vertexNumber+1):
            for j in range(len(tableOfEdges[0])):
                if tableOfEdges[j][i] != 0:
                    L = L+tableOfEdges[j][i].edgelength
        return L

    def calculateVertexCurvature(self,vertex,tableOfEdges):
        vertexCurvature = 0
        for i in range(len(tableOfEdges[0])):
            if tableOfEdges[vertex][i] != 0:
                vertexCurvature = vertexCurvature+tableOfEdges[vertex][i].edgecurvature
            elif tableOfEdges[i][vertex] != 0:
                vertexCurvature = vertexCurvature+tableOfEdges[i][vertex].edgecurvature
        vertexCurvature = vertexCurvature/2
        return vertexCurvature


    def checkLCSC(self,tableOfEdges,error=.0001):
        LCSC = True
        listOfLengths = []
        for i in range(self.vertexNumber+1):
            listOfLengths.append([0])
        for i in range(len(self.edgeList)):
            listOfLengths[self.edgeList[i][0]].append(self.edgetable[self.edgeList[i][0]][self.edgeList[i][1]].edgelength)
            listOfLengths[self.edgeList[i][1]].append(self.edgetable[self.edgeList[i][0]][self.edgeList[i][1]].edgelength)
       # for i in range(self.vertexNumber):
       #     print(i,self.calculateVertexCurvature(i,tableOfEdges),self.LEHR * L)
        for i in range(1,self.vertexNumber+1):
            if math.fabs(self.calculateVertexCurvature(i,tableOfEdges)-(self.LEHR * sum(listOfLengths[i])/2)) > (error):
                #print(math.fabs(self.calculateVertexCurvature(i,tableOfEdges)-(self.LEHR * L)))
                LCSC = False
        return LCSC

    def checkLEinstein(self,error=.0001):
        LEinstein = True
        for i in range(len(self.edgeList)):
            curvature=self.edgetable[self.edgeList[i][0]][self.edgeList[i][1]].edgecurvature
            LEHRl=self.edgetable[self.edgeList[i][0]][self.edgeList[i][1]].edgelength
            LEHRl=LEHRl*self.LEHR
            if math.fabs(curvature-LEHRl)> error:
                LEinstein= False
        return LEinstein

    #inpout: table of edges and number of verticies
    #output: none, creates a list such that the n-1 item of the list corresponds with the sum of all the edges incident on edge n
    #Author: EM, 02/04/15
    #change log:

    def getEdgeSums(self,tableOfEdges,numberOfVertices):
        edgeSums = []
        for i in range(1,numberOfVertices+1):
            sums = 0
            for j in range(1,numberOfVertices+1):
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
    def __init__(self,conformalVariations, backgroundFile = 'backgroundMetric.txt', manifoldFile = 'manifoldExample.txt'):
        self.edgetable = []#A list of lists of edges such that edge x,y is stored in edgetable[x][y]
        self.tetrahedralist = []#A list of tetrahedron objects
        self.edgeList = []#A list of edge names where edgeList[i]=[[a],[b]]
        self.vertexNumber = 0
        self.LEHR = 10000
        self.LCSC = False
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
        self.vertexNumber = self.createTetrahedraList(tetrahedron)
        self.createEdgeTable(self.edgetable,self.vertexNumber)
        self.fillEdgeTable(self.tetrahedralist,self.edgetable,self.edgeList,backgroundFile,conformalVariations)
        self.sumOfEdgesAtVertex = self.getEdgeSums(self.edgetable,self.vertexNumber)
        illegalTetrahedrons = self.advancedCheckLegalTetrahedra(self.tetrahedralist,self.edgetable)
        if illegalTetrahedrons == True:
            self.good = False
        else:
            for i in range(len(self.tetrahedralist)):
                self.tetrahedralist[i].initalizeTriangles(self.edgetable)
            self.totalLength = self.findTotalEdgeLength(self.edgetable)
            for i in range(len(self.tetrahedralist)):
                self.tetrahedralist[i].calculateDihedralAngles(self.edgetable)
                self.tetrahedralist[i].getTetCenDis()
            for i in range(len(self.edgeList)):
                self.edgetable[self.edgeList[i][0]][self.edgeList[i][1]].calculateEdgeCurvature(self.tetrahedralist)
                self.edgetable[self.edgeList[i][0]][self.edgeList[i][1]].calculateEdgeLengthStar(self.tetrahedralist)
            self.edgeStarOverEdgeTotal = self.getEdgeStarOverEdgeSums(self.edgetable,self.vertexNumber)
            self.LEHR = self.findLEHR(self.edgeList,self.edgetable)
            #self.showEdgeTable(self.edgetable)
            for i in range(self.vertexNumber):
                self.vertexCurvatureList.append(self.calculateVertexCurvature(i+1,self.edgetable))
            self.LCSC = self.checkLCSC(self.edgetable)
            self.LEinstein = self.checkLEinstein()
            self.good = True
    # 10/16/14 Michael; replaced for loop to check legal tetrahedron to advancedcheckLegalTetrahedron command


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


def hessianSame(metric,vertex):
    partial = 2*metric.edgeStarOverEdgeTotal[vertex-1]/metric.totalLength
    partial = partial + metric.vertexCurvatureList[vertex-1]/(2*metric.totalLength)
    partial = partial - metric.LEHR*metric.sumOfEdgesAtVertex[vertex-1]**2/(4*metric.totalLength)
    partial = partial - metric.sumOfEdgesAtVertex[vertex-1]*metric.vertexCurvatureList[vertex-1]/(metric.totalLength**2)
    partial = partial + metric.LEHR*metric.sumOfEdgesAtVertex[vertex-1]**2/(2*metric.totalLength**2)
    return partial
def hessianDifferent(metric,vertex1,vertex2):
    if vertex1 > vertex2:
        lijStarOverlij = metric.edgetable[vertex2][vertex1].edgelengthStar/metric.edgetable[vertex2][vertex1].edgelength
        edgeCurvature = metric.edgetable[vertex2][vertex1].edgecurvature
        edgeLength = metric.edgetable[vertex2][vertex1].edgelength
    else:
        lijStarOverlij = metric.edgetable[vertex1][vertex2].edgelengthStar/metric.edgetable[vertex1][vertex2].edgelength
        edgeCurvature = metric.edgetable[vertex1][vertex2].edgecurvature
        edgeLength = metric.edgetable[vertex1][vertex2].edgelength
    partial = -2*lijStarOverlij/metric.totalLength
    partial = partial + edgeCurvature/(4*metric.totalLength)
    partial = partial - metric.LEHR*edgeLength/(4*metric.totalLength)
    partial = partial - metric.vertexCurvatureList[vertex2-1]*metric.sumOfEdgesAtVertex[vertex1-1]/(2*metric.totalLength**2)
    partial = partial + metric.sumOfEdgesAtVertex[vertex1-1]*metric.LEHR*metric.sumOfEdgesAtVertex[vertex2-1]/(2*metric.totalLength**2)
    partial = partial - metric.sumOfEdgesAtVertex[vertex2-1]*metric.vertexCurvatureList[vertex1-1]/(2*metric.totalLength**2)
    return partial

def makePosDef(target,add = 1):
    eigns = numpy.linalg.eigvals((target))
    smallEig = min(eigns)
    if(smallEig <= 0):
        target = target + (math.fabs(smallEig)+add)*numpy.identity(target.size**.5)
    eigns2=numpy.linalg.eigvals(target)
    return target

def newtonsMethod(mainMetric,convar ,background,triagulation,stepSize = 1,gradPos = True,numberCalls = 0):
    gamma = 0
    grad = []
    hessian = []
    if(mainMetric.good == False):
        return -1
    else:
        for i in range(mainMetric.vertexNumber):
            #creating the gradient
            temp = mainMetric.vertexCurvatureList[i]
            temp = temp-mainMetric.LEHR*.5*mainMetric.sumOfEdgesAtVertex[i]
            temp = temp/mainMetric.totalLength
            grad.append(-temp)
        for i in range(mainMetric.vertexNumber):
            #creating the hessian
            hessian.append([0]*mainMetric.vertexNumber)
            for j in range(mainMetric.vertexNumber):
                if mainMetric.edgetable[i+1][j+1] == 0 and mainMetric.edgetable[j+1][i+1] == 0 and i!=j:
                    if mainMetric.edgetable[i+1][j+1] == 0:
                        hessian[i][j] = 0
                elif i == j:
                    hessian[i][j] = hessianSame(mainMetric,i+1)
                else:
                    hessian[i][j] = hessianDifferent(mainMetric,i+1,j+1)
        question = hessian
        questionTemp = numpy.array(question)
        newquestion = makePosDef(questionTemp,1)
        question = newquestion.tolist()
        #tempQuestion = copy.deepcopy(question)
        checkPosDef  = numpy.array(question)
        if(min(numpy.linalg.eigvals(checkPosDef)) <= 0):
            print("NOT POS DEF!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!,")
            question=makePosDef(newquestion).tolist()
        #Make hessian pos def
        for i in range(len(grad)):
             question[i].append(grad[i])
        ToReducedRowEchelonForm(question)
        answer = [0]*len(grad)
        temp = []
        for i in range(len(grad)):
            answer[i] = stepSize*question[i][len(question[i])-1]
        #while(i != i):
        newConvar = []
        for i in range(len(convar)):
            newConvar.append(convar[i])
            newConvar[i] = newConvar[i]+answer[i]
        temp = metric(newConvar, background, triagulation)
        transposeTimesSk = 0
        counter = 0
        for i in range(len(convar)):
            transposeTimesSk = transposeTimesSk + convar[i]*grad[i]
        while(temp.LEHR >= mainMetric.LEHR+gamma*transposeTimesSk and stepSize > 1/2**10 and temp.LCSC == False):
            stepSize = stepSize/2
            counter = counter+1
            #if(stepSize < 1/2**14):
            for i in range(len(convar)):
                #print(stepSize)
                newConvar[i] = convar[i]
                newConvar[i] = newConvar[i]+answer[i]*stepSize
            temp = metric(newConvar,background,triagulation)
        #print(mainMetric.LEHR-temp.LEHR)

        return newConvar

def modifyBackground(c1,c2,filename):
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
    lengthList = [(math.exp(c1))**.5,(math.exp(c2))**.5,(math.exp(-c1-c2))**.5,(math.exp(-c1-c2))**.5,(math.exp(c2))**.5,(math.exp(c1))**.5]
    backgroundMetric = open(filename,"w")
    for i in range(len(nameList)):
        backgroundMetric.write(str(nameList[i][0])+","+str(nameList[i][1])+"\n"+str(lengthList[i])+"\n")
    backgroundMetric.close()

def doubleTetrahedronWalk(numberVertices,backgroundfile,triangulation,restarts = 100,numberBackgrounds = 10):
    results = []
    failures = []
    working = True
    results.append(["c1","c2","f1","f2","f3","f4","LEHR","LCSC","L-Einstein","numberRestarts"])
    failures.append(["c1","c2","Conformal Variations"])
    for k in range(numberBackgrounds):
        c1 = random.random()*10-5
        c2 = random.random()*10-5
        print("background number "+str(k+1) +" out of " + str(numberBackgrounds))
        while(not (-math.exp(c2)+math.exp(c1)+math.exp(-c1-c2) > 0 and math.exp(c2)+math.exp(c1)-math.exp(-c1-c2) > 0 and -math.exp(c2)+math.exp(c1)-math.exp(-c1-c2) < 0)):
            c1 = random.random()*10-5
            c2 = random.random()*10-5
        for j in range(restarts):
            working = True
            happyConVar = False
            while(happyConVar == False):
                conVar = []
                modifyBackground(c1,c2, "backgroundMetric.txt")
                for i in range(numberVertices):
                    #sets conformal variations to 0
                    #conVar.append(0)
                    #sets random conformal variations
                    conVar.append(random.random())
                orignalConVar = []
                orignalConVar = copy.deepcopy((conVar))
                test = metric(conVar,backgroundfile,triangulation)
                happyConVar = test.good
            if(test.good == False):
                    print("Illegal Initial")
            ConformalStep = newtonsMethod(test,conVar,backgroundfile,triangulation)
            lastLEHR = test.LEHR
            conVarStore = conVar
            stepSize = 2
            newConformal = []
            gradPos = True
            conVarStore = conVar
            for i in range(5000):
                newConformal = newtonsMethod(test,conVar,backgroundfile,triangulation,1,gradPos)
                if(newConformal == -1):
                    working = False
                    break
                test = metric(newConformal,backgroundfile,triangulation)
                ConVarStore = copy.deepcopy(conVar)
                conVar = newConformal
                if i%1000 == 0 and i != 0:
                    print(test.LEHR)
                    print(i)
                if test.LCSC == True:
                    print("Found at step "+str(i))
                    working = True
                    foundAt = j
                    break
            if working == True:
                results.append([c1,c2,conVar[0],conVar[1],conVar[2],conVar[3],test.LEHR,test.LCSC,test.LEinstein,foundAt])
                print("Found at restart "+str(j))
                break
        if working == False:
            failures.append([c1,c2])
    return [results,failures]


def main():
   storage = str(0)+".txt"
   steps = 5000
   restarts = 1
   LEHRList = []
   numberVertices=4
   Restarts =25 #number of new sets of conformal variations tested
   numberOfBackgrounds=2
   #seed=4741252
   #seed=263594
   seed=56932684
   #seed=9865721
   triangulation='manifoldExample3.txt'
   backgroundfile='backgroundMetric.txt'
   faceInfo = " "
   print("Hello World!\n")
   random.seed(seed)
   all = doubleTetrahedronWalk(numberVertices,backgroundfile,triangulation,Restarts,numberOfBackgrounds)
   notFounds = all[1]
   store = all[0]
   for i in range(len(store)):
       print(store[i])
   print("\n\n")
   for i in range(len(notFounds)):
       (notFounds[i])
   resultsFile = open("Results.txt","w")
   for i in range(len(store)):
       for j in range(len(store[i])):
            resultsFile.write(str(store[i][j])+ "  ")
       resultsFile.write("\n")
   for i in range(len(notFounds)):
        for j in range(len(notFounds[i])):
            resultsFile.write(str(notFounds[i][j])+ "  ")
        resultsFile.write("\n")
   resultsFile.close()


main()
