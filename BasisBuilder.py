#10/16/14 Michael: Added advancedCheckLegalTetrahedra function for checking legal tetrahedron and used it main, made minor change checkLegalTetrahedron
import math
import numpy as np
import copy

edgetable = []#A list of lists of edges such that edge x,y is stored in edgetable[x][y]
tetrahedralist = []#A list of tetrahedron objects
edgeList = []#A list of edge names where edgeList[i]=[[a],[b]]
vertexNumber = 0



class Edge:
    #vertex1: an integer that corresponds with the first part of the name of an edge
    #vertex2: an integer that corresponds with the second part of the name of an edge
    #edgelength: a number that represents how long the edge is
    #edgecurvature: a number that is edge curvature
    #tetrahdraEdgeIsIn: location in tetrahedraList of all tetrahedra that contain that edge

    #__init__(self,vertex1 = 1,vertex2 = 2): returns nothing, sets default values
    #calculateEdgeCurvature(self,listOfTetrahedra): returns nothing, sets the edge curvature

    #input: two numbers that represent the name of an edge
    #output: none, sets default values of members of the edge class
    #Author: Prof. Young, 10/2/2014
    #change log: Erin, Lincoln, Michael, Tyler 10/3/2014
    def __init__(self,vertex1 = 1,vertex2 = 2):
        self.vertex1 = vertex1
        self.vertex2 = vertex2
        self.edgelength = 1
        self.edgecurvature = 0
        self.tetrahedraEdgeIsIn = []
    #10/3/2014: added tetrahedraEdgeIsIn

    #input: list of Tetrahedron objects
    #output: none, sets edgecurvature to the curvature of the edge
    #author: MELT, 10/14/2014
    #change log: none
    def calculateEdgeCurvature(self,listOfTetrahedra):
        diList = []
        for i in range (len(self.tetrahedraEdgeIsIn)):
            tetLocation = self.tetrahedraEdgeIsIn[i]
            # Finds the dihedral angle of edge in tetrahedron at tetrahedraEdgeIsIn[i]
            singleDi = listOfTetrahedra[tetLocation].dihedralanglelist[listOfTetrahedra[tetLocation].edgesintetrahedron.index([self.vertex1,self.vertex2])]
            diList.append(singleDi)
        self.edgecurvature = (2*math.pi-sum(diList))*self.edgelength # Dihedral angle formula
    #
        
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
        self.edgesintetrahedron = []  #stored as lists ex, edge 1,2 stored as [1,2]
        self.dihedralanglelist = []
    #

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
        if det<0 or det==0:
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

    #input: Eij,Eik,Eil,Ejk,Ejl,Ekl, numbers that represent edge lengths of a tetrahedron
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
        self.dihedralanglelist.append(math.acos(self.calDiAngle(edge14,edge12,edge13,edge24,edge34,edge23)))
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


#input: number of vertices, currently hard coded at 15
#output: none, makes the edge table of proper size, i.e. 15+1
#author: MELT, 10/3/14
#change log:
def createEdgeTable(numberOfVertices = 5):
    numberOfVertices = numberOfVertices+1
    for row in range(numberOfVertices):
        edgetable.append([])
        for column in range(numberOfVertices):
            edgetable[row].append(0)
#


#input: table of edge objects
#output: none, debugging function
#author: MELT, 10/7/2014
#change log:
def showEdgeTable(tableOfEdges):
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
def createTetrahedraList(primalList):
    numberOfVertices = 0
    for i in range(len(primalList)):
        tetrahedralist.append(Tetrahedron(primalList[i][0],primalList[i][1],primalList[i][2],primalList[i][3]))
        tetrahedralist[i].edgesintetrahedron.append([primalList[i][0],primalList[i][1]])
        tetrahedralist[i].edgesintetrahedron.append([primalList[i][0],primalList[i][2]])
        tetrahedralist[i].edgesintetrahedron.append([primalList[i][0],primalList[i][3]])
        tetrahedralist[i].edgesintetrahedron.append([primalList[i][1],primalList[i][2]])
        tetrahedralist[i].edgesintetrahedron.append([primalList[i][1],primalList[i][3]])
        tetrahedralist[i].edgesintetrahedron.append([primalList[i][2],primalList[i][3]])
        if numberOfVertices < tetrahedralist[i].vertex4:
            numberOfVertices = tetrahedralist[i].vertex4
    return numberOfVertices
#10/16/2014 MELT added vertex counter

#

def findSpanTree(spanM,size):
    numberEdges = sum(sum(spanM))/2
    sT = np.array([[0]*size]*size)
    found = [False]*size
    explore = []
    for i in range(len(spanM)):
        if spanM[0][i] == 1:
            explore.insert(0,[0,i])
    found[0] = True
    while len(explore) != 0:
        focus = explore.pop()
        if found[focus[1]] == False:
            sT[focus[0]][focus[1]] = 1
            sT[focus[1]][focus[0]] = 1
            focus = focus[1]
            found[focus] = True
            for i in range(len(spanM)):
                if spanM[focus][i] == 1:
                    explore.insert(0,[focus,i])
    print(sT)
    return sT

def createOddCyle(sT,graph):
    oddCycle = False
    legalEdge = False
    size = len(graph)
    count = 0
    while oddCycle == False:
        while legalEdge == False:
            j = count%size
            i = math.floor(count/size)
            if sT[i][j] == 0 and graph[i][j] == 1:
                legalEdge = True
            count = count+1
        print("Add edge")
        sT[i][j] = 1
        sT[j][i] = 1
        print("Check cycles")
        oddCycle = findCycleOdd(sT,i,j)
        if(oddCycle == False):
            sT[i][j] = 0
            sT[j][i] = 0
    return sT

def findCycleOdd(tree,start,end):
    print("Build Storage")
    found = False
    edgesUsed = [[0 for i in range(len(tree))] for j in range(len(tree))]
    edgesUsed[start][end] = 1
    edgesUsed[end][start] = 1
    print("Start Search")
    cyle = searchForCycle(start,end,tree,edgesUsed,False,True)
    print("check if cycle is odd")
    if len(cyle)%2 == 1:
        found = True
    print("Return if found or not")
    return found

def searchForCycle(start,end,tree,explored,wantEven,startEven = True):
    cycleList = [-1]*len(tree)
    cycleList[start] = end
    found = DFS(start,tree,explored,end,cycleList,wantEven,startEven)
    if found == True:
        cycle = []
        cycle.insert(0,end)
        cycle.insert(0,cycleList[end])
        while(cycle[0] != start and len(cycle) <= len(tree)):
            cycle.insert(0,cycleList[cycle[0]])
        if cycle[0] != start:
            cycle = [-1]
            found = False
    else:
        cycle = [-1]
    if found == False:
        print("No cycle present")
    return cycle

def DFS(start,tree,explored,end,cycleList=[],wantEven = True,even = False):
    good = False
    for i in range(len(tree[start])):
        if tree[start][i] == 1 and explored[start][i] != 1 and i == end and even == wantEven:
            cycleList[i] = start
            good = True
            break
        elif tree[start][i] == 1 and explored[start][i] != 1 and i != end:
            cycleList[i] = start
            explored[start][i] = 1
            explored[i][start] = 1
            good = DFS(i,tree,explored,end ,cycleList,wantEven,not even)
            explored[start][i] = 0
            explored[i][start] = 0
            if good == True:
                break
    return good
#input: none
#output: none, super-mega-function that does EVERYTHING! prints LEHR
#author: METAL, 10/8/2014
#change log: Michael 10/16/4
def main():
    print("Hello World")
    readFile = open('manifoldExample4.txt')
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
    print(tetrahedron)
    vertexNumber = createTetrahedraList(tetrahedron)
    adjMatrix = np.array([[0]*vertexNumber]*vertexNumber)
    print("Create graph")
    for i in range(len(tetrahedron)):
        vert1 = tetrahedron[i][0]-1
        vert2 = tetrahedron[i][1]-1
        vert3 = tetrahedron[i][2]-1
        vert4 = tetrahedron[i][3]-1
        adjMatrix[vert1][vert2]=1
        adjMatrix[vert1][vert3]=1
        adjMatrix[vert1][vert4]=1
        adjMatrix[vert2][vert3]=1
        adjMatrix[vert2][vert4]=1
        adjMatrix[vert3][vert4]=1

        adjMatrix[vert2][vert1]=1
        adjMatrix[vert3][vert1]=1
        adjMatrix[vert4][vert1]=1
        adjMatrix[vert3][vert2]=1
        adjMatrix[vert4][vert2]=1
        adjMatrix[vert4][vert3]=1
    print(adjMatrix)
    print("Create spanning tree")
    spanTree = findSpanTree(adjMatrix,vertexNumber)
    print("Add edge to create odd cycle")
    oddCycle = createOddCyle(spanTree,adjMatrix)
    print(oddCycle)
    print("look at dendroid")
    dendroid = np.subtract(adjMatrix,oddCycle)
    print(dendroid)
    #repeate for each edge in dendroid
    for i in range(len(dendroid)):
        for j in range(len(dendroid)):
            if dendroid[i][j] == 1 and j > i:
                oddCycle[i][j] = 1
                oddCycle[j][i] = 1
                #awesome math stuff
                explored = [[0 for i in range(len(spanTree))] for j in range(len(spanTree))]
                explored[i][j] = 1
                explored[j][i] = 1
                cycle = searchForCycle(i,j,oddCycle,explored,True,True)
                print(str(i)+str(cycle)+str(j))
                oddCycle[i][j] = 0
                oddCycle[j][i] = 0
    print("find even cycle/two odd cycles")
    print("assign values")


    createEdgeTable(vertexNumber)


main()
