#10/16/14 Michael: Added advancedCheckLegalTetrahedra function for checking legal tetrahedron and used it main, made minor change checkLegalTetrahedron
import math

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
def createEdgeTable(numberOfVertices = 15):
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


# Input: list of tetrahedra and empty table of edges and a name of file to be created, default backgroundMetric
# Output: fills the empty table of edges and outputs a file of edge names in order they were found with edge lengths associated with the edge directly below
# Author, M, 10/27/14
# Change Log
def getBackgroundManifoldFile(listOfTetrahedra,tableOfEdges,fileName = "backgroundMetric.txt"):
    outfile = open(fileName,"w")
    for i in range(len(listOfTetrahedra)):
        for j in range(len(listOfTetrahedra[i].edgesintetrahedron)):
            # If edge is not in the the tableOfEdges, add it
            if tableOfEdges[listOfTetrahedra[i].edgesintetrahedron[j][0]][listOfTetrahedra[i].edgesintetrahedron[j][1]] == 0:
                # Assigns edge its name
                tableOfEdges[listOfTetrahedra[i].edgesintetrahedron[j][0]][listOfTetrahedra[i].edgesintetrahedron[j][1]] = Edge(listOfTetrahedra[i].edgesintetrahedron[j][0],listOfTetrahedra[i].edgesintetrahedron[j][1])
                # Adds tetrahedran to list of tetrahedra edge is in
                tableOfEdges[listOfTetrahedra[i].edgesintetrahedron[j][0]][listOfTetrahedra[i].edgesintetrahedron[j][1]].tetrahedraEdgeIsIn.append(i)
                outfile.write(str(listOfTetrahedra[i].edgesintetrahedron[j][0])+","+str(listOfTetrahedra[i].edgesintetrahedron[j][1])+"\n1\n")
            # adds edge to tetrahera list if edge already exists
            else:
                tableOfEdges[listOfTetrahedra[i].edgesintetrahedron[j][0]][listOfTetrahedra[i].edgesintetrahedron[j][1]].tetrahedraEdgeIsIn.append(i)
    outfile.close()
#


#input: none
#output: none, super-mega-function that does EVERYTHING! prints LEHR
#author: METAL, 10/8/2014
#change log: Michael 10/16/4
def main():
    print("Hello World")
    readFile = open('manifoldExample3.txt')
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
    vertexNumber = createTetrahedraList(tetrahedron)
    createEdgeTable(vertexNumber)
    getBackgroundManifoldFile(tetrahedralist,edgetable)


main()
