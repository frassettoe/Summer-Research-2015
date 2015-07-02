__author__ = 'Michael'


import numpy as np
import random as rand
import copy
import warnings as warning
import math


#In this code, all graphs are assumed to be in matrix form and stored as numpy arrays
#Convention: edge i,j is such that i <= j
#When input is graph, that means send the full graph (not dendroid/spanning tree)

#Input: Size n by n graph, goal of evantually creating a random graph with percentFull percent of maximum edge lengths
#Output: graph
#Author: M, 6/30/15
#Change Log: none
def buildGraph(size,percentFull = 1):  #percentful not used...yet
    print("Build Graph")
    graph = np.zeros([size,size])  #creates empty graph
    print(graph)
    for i in range(size):  #adds every edge possible
        for j in range(size):
            if i != j:
                graph[i,j] = 1
    print(graph)
    return graph

#Input: graph, in matrix form
#Output: Spanning tree of graph
#Author: M, 6/30/15
#Changes: none
def buildSpanningTree(graph):
    print("Build Spanning Graph")
    spanTree = np.zeros([len(graph),len(graph)])  #creates empty tree, same size as graph
    for i in range(len(graph)):  #builds spanning tree
        for j in range(len(graph[i])):
            if graph[i][j] == 1 and i <= j: #if edge exists
                if len(findSingleCyclePath(spanTree,i,j)) == 1:  #if ending of edge is not in the graph
                    spanTree[i][j] = 1  #  add edge to graph
                    spanTree[j][i] = 1
    return spanTree
# Early attempt at finding all cycles in a graph, never called, not tested
# Author: M, 6/30/15
# def findCycles(graph):  #never called
#     print("Find fundimental Cycles")
#     spanTree = buildSpanningTree(graph)
#     cycleMatrix = []
#     print(spanTree)
#     chords = np.subtract(graph,spanTree)
#     print("Chords")
#     print(chords)
#     print("Cycles")
#     for i in range(len(graph)):
#         for j in range(len(graph)):
#             if chords[i,j] == 1:
#                 spanTree[i,j] = 1
#                 spanTree[j,i] = 1
#                 chords[i,j] = 0
#                 chords[j,i] = 0
#                 cycle = findSinglePath(spanTree,i,j)
#                 spanTree[i,j] = 0
#                 spanTree[j,i] = 0
#                 cycleMatrix.append([0]*len(graph))
#                 for i in range(len(cycle)):
#                     cycleMatrix[-1][cycle[i]] = 1
#     cycleMatrix = np.array(cycleMatrix)
#     return cycleMatrix

#Input: graph, starting node, ending node
#output: List of nodes forming a path between start and end
#Author: M, 06/30/15
#Change log: need to check what occurs when no path exists
def findSingleCyclePath(graph,start,end):
    reAdd = False
    parents = [-1]*len(graph) #list of parents of nodes, aka parents[3] is the node that connects to node 3 and -1 represents node not visited
    parents[start] = end #Notes that we waant a loop, and in the limits of the graphs we use, start and end are almost always connected, when they are not connected, still correct path is found
    if graph[start,end] == 1:
        graph[start,end] = 0 #removes connection between start and end (otherwise the path found would just be start end
        graph[end,start] = 0
        reAdd = True
    toLookAt = [] #stack of edges to visit
    for i in range(len(graph)): #inital fill of toLookAt, look at all edges connected to node start
        if graph[start,i] == 1: #if edge exists
            toLookAt.append([start,i]) #add it to toLookAt
    while parents[end] == -1 and len(toLookAt) != 0: #while ending node not found and toLookAt is not empty
        focus = toLookAt.pop() #remove last item added to toLookAt
        if parents[focus[1]] == -1: #if node visiting has not already been visited
            parents[focus[1]] = focus[0] #updates paraent for the node just visited
            for i in range(len(graph)): #look at each possible edge connected to node
                if graph[focus[1],i] == 1: #if edge exists
                    toLookAt.append([focus[1],i]) #add edge to toLookAt
    cycle = [end] #adds end of path to cycle
    if reAdd == True:
        graph[start,end] = 1 #adds edge removed earlier to graph
        graph[end,start] = 1
    if parents[end] != -1:
        while cycle[-1] != start and len(cycle) < len(graph)+1: #while the start of the path has not been found and the length of the cycle found is not longer than the number of nodes (this done to avoid infinite loop)
            cycle.append(parents[cycle[-1]])  #add parent of last edge in cycle
        if len(cycle) >= len(graph)+1: #warn due to untested case
            warning.warn("No path Found, possible fatal error")
    return cycle
#Input: graph, dendroid, spanning tree
#Output: odd cycle in spanning tree with addition of edge from graph not in present spanning tree
#Author: M, 6/30/15
#Change Log: none
def findOddCyclesFromSpanningTree(graph,dendroid,sT):
    cycles = []
    for i in range(len(graph)): #for each possible edge in graph
        for j in range(len(graph)):
            if i < j and dendroid[i,j] == 1:  #if edge follows convention and is in dendroid
                sT[i,j] = 1 #add edge to spanning tree
                sT[j,i] = 1
                cycles = findSingleCyclePath(sT,i,j) #find fundemental cycle of spanning tree and added edge
                addedEdge = [i,j]
                sT[i,j] = 0 #removes edge added earilier from spanning tree
                sT[j,i] = 0
                if len(cycles)%2 == 1:  #if cycle is odd
                    break #done
        if len(cycles)%2 == 1: #if cycle is odd
            break #done, condition check doubled due to nested for loop
    dendroid[addedEdge[0],addedEdge[1]] = 0 #removes edge used to create odd cycle from dendroid (making dendroid an actual dendroid)
    dendroid[addedEdge[1],addedEdge[0]] = 0
    print("Dendroid is:")
    print(dendroid)
    print("Extra Edge is")
    print(addedEdge)
    return cycles
#input: graph, dendroid, spanning tree, cycles (style two: primal fundemental cycle position 0 and mutated fundemental cycle position 1)
#output: Basis of graph
#Author: M, 06/30/15
#Changes: none
def buildBasisWithEvenCycles(graph,dendroid,sT,cycles2):
    vectors = []
    for i in range(len(graph)):  #for every possible edge in graph
        for j in range(len(graph)):
            if i < j and dendroid[i,j] == 1:  #if edge follows convention and is in dendroid
                sT[i,j] = 1  #add edge from dendroid to spanning tree
                sT[j,i] = 1
                cycles2[1] = findSingleCyclePath(sT,i,j) #find fundemental cycle (called mutated fundemental cycle)
                sT[i,j] = 0  #removes edge from dendroid from spanning tree
                sT[j,i] = 0
                print(str(i)+" , " +str(j))
                vectors.append(getVector(cycles2,graph,i,j)) #turns cycle into vector and adds it to vector list
    return vectors

#Input: cycle (starting with i of dendroid edge and ending with j of dendroid edge) and graph
#output: vector from cycle and graph
#Author: ME, 06/31/15
#Change log: none
def getCycleVector(cycle,graph):
    vector = []
    cycleGraph = np.zeros([len(graph),len(graph)]) #create graph that just contains edges of cycles
    for i in range(len(cycle)):  #for each node in cycle
        if cycle[i] < cycle[i-1]: #if convention is followed by adjacent nodes of cycle
            cycleGraph[cycle[i],cycle[i-1]] = 1-2*(i%2) #add node to cycle graph (switching between adding a 1 and -1)
        else: #otherwise due same thing with inputs switched so convention followed
            cycleGraph[cycle[i-1],cycle[i]] = 1-2*(i%2)
    for i in range(len(graph)):  #for every edge in graph, this is done to ensure edge i,j appears in the same position for every vector of the same graph
        for j in range(len(graph[i])):
            if i < j and cycleGraph[i,j] != 0: #if convention followed and edge in cycle
                vector.append(int(cycleGraph[i,j])) #add one or zero to vector
            elif i < j and graph[i,j] != 0: #else if convention followed and edge in graph
                vector.append(0) #add zero
    return vector
#Input: primal cycle, mutated cycle, graph, node that is in both cycles
#output: vector
#Author: ME 06/31/15
#Change Log: none
def getCycleCycleVector(cycle1,cycle2,graph,overlap):
    vector = []
    cycle = []
    pointOfOverlapInCycle1 = cycle1.index(overlap)  #find the location of overlap in the primal cycle
    #combine cycle 1 and cycle 2
    for i in range(len(cycle2)): #for each node in mutated cycle
        if cycle2[i] == overlap: #if overlap occurs
            for j in range(len(cycle1)): #add nodes from primal cycle
                cycle.append(cycle1[(pointOfOverlapInCycle1+j)%len(cycle1)])
        cycle.append(cycle2[i]) #add node from mutated cycle
    cycleGraph = np.zeros([len(graph),len(graph)])
    #add cycle to cycle graph, switching between 1 and -1
    for i in range(len(cycle)):
        if cycle[i] < cycle[i-1]:
            cycleGraph[cycle[i],cycle[i-1]] = 1-2*(i%2)
        else:
            cycleGraph[cycle[i-1],cycle[i]] = 1-2*(i%2)
    #creates vector, following pattern of getCycleVector
    for m in range(len(graph)):
        for n in range(len(graph[m])):
            if m < n and cycleGraph[m,n] != 0:
                vector.append(int(cycleGraph[m,n]))
            elif m <  n and graph[m,n] != 0:
                vector.append(0)
    return vector

def getCyclePathCycle(cycle1,cycle2,path,graph):
    warning.warn("Untested Function Called")
    vector = []
    cycle = []
    value = 1
    pointOfOverlapInCycle1WithPath = cycle1.index(path[0])  #find the location of overlap in the primal cycle
    pointOfOverlapInCycle2WithPath = cycle2.index(path[-1])  #find the location of overlap in the primal cycle
    #combine cycle 1 and cycle 2
    for i in range(len(cycle2)): #for each node in mutated cycle
        if i == pointOfOverlapInCycle2WithPath: #if overlap occurs
            for k in range(len(path)):
                cycle.append(path[len(path)-k-1])
            for j in range(len(cycle1)-1): #add nodes from primal cycle
                cycle.append(cycle1[(pointOfOverlapInCycle1WithPath+j+1)%len(cycle1)])
            for k in range(len(path)-1):
                cycle.append(path[k])
        cycle.append(cycle2[i]) #add node from mutated cycle
    cycleGraph = np.zeros([len(graph),len(graph)])
    #add cycle to cycle graph, switching between 1 and -1
    for i in range(len(cycle)):
        if cycle[i] < cycle[i-1]:
            cycleGraph[cycle[i],cycle[i-1]] += value
        else:
            cycleGraph[cycle[i-1],cycle[i]] += value
        value = -value
    #creates vector, following pattern of getCycleVector
    for m in range(len(graph)):
        for n in range(len(graph[m])):
            if m < n and cycleGraph[m,n] != 0:
                vector.append(int(cycleGraph[m,n]))
            elif m <  n and graph[m,n] != 0:
                vector.append(0)

    return vector

#Input: cycle, graph, two nodes from dendroid (checks if mutated cycle is even or odd)
#Output: vector
#Author: ME, 06/31/15
#Change Log: Updated terminology of prints and removed  line vector = [0]*len(graph) from start of graph (M, 07/1/15)
def getVector(cycles,graph,nodeFromDendroid1,nodeFromDendroid2):
    if len(cycles[1])%2 != 0:  #mutated fundemental not cyle even
        print("Mutated funedemntal cycle odd")
        vector = newCycleOddSoFindEven(graph,cycles,nodeFromDendroid1,nodeFromDendroid2) #call fucntion to consider this
    else:  #Mutated fundemental cycle even, yay!
        print("Mutated Create Vector, new fundemental cycle even")
        print(cycles)
        vector = getCycleVector(cycles,graph)
        #1. turn cycle into graph  #notes on plans for turning cycle into vector!
        #2. use standard means of going from graph to list
        # for i in range(len(cycles)):
        #     vector[cycles[1][i]] = 1-2*(i%2)
    return vector

#input: graph, cycles style 2 (both odd), start and end of desired cycle
#ouput: vector, following rules defined to us
#Author: ME 06/31/15
#Change Log: none
def newCycleOddSoFindEven(graph,cycles,i,j):
    edgeTable = np.zeros([len(graph),len(graph)])  #graph that contains both cycles
    #add both fundemental cycles
    for k in range(len(cycles[0])): #adds primal cycle
        edgeTable[cycles[0][k],cycles[0][k-1]] = 1
        edgeTable[cycles[0][k-1],cycles[0][k]] = 1
    for k in range(len(cycles[1])): #adds mutated cycle
        if edgeTable[cycles[1][k],cycles[1][k-1]] == 1:  #remove edges that appear in both cycles
            edgeTable[cycles[1][k],cycles[1][k-1]] = 0
            edgeTable[cycles[1][k-1],cycles[1][k]] = 0
        else:
            edgeTable[cycles[1][k],cycles[1][k-1]] = 1   #otherwise add edge
            edgeTable[cycles[1][k-1],cycles[1][k]] = 1
    #print(edgeTable)
    cycles3 = findSingleCyclePath(edgeTable,i,j)  #Find cycle from two dendroids (not sure what will happen if cycles are not connected)
    print(cycles3)
    if len(cycles3)%2 != 0:  #if new cycle found is not even
        print("New combined cycle Odd")
        vector = advancedOddSearch(graph,cycles) #check one of three possible cases
    else:
        vector = getCycleVector(cycles3,graph) #otherwise we have even cycle!
        print("Combined cycle even")
    return vector

def advancedOddSearch(graph,cycles):
    #check to see if cycles overlap on nodes once or multiple times
    overlap = -1
    overlapTwiceOrMore = False
    for m in range(len(cycles)):
        for n in range(len(cycles[m])):
            if overlap != -1 and cycles[0][m] == cycles[1][n]:
                overlapTwiceOrMore = True
                break
            if cycles[0][m] == cycles[1][n]:
                overlap = cycles[0][m] #stores point of overlap
        if overlapTwiceOrMore == True:
            break
    if overlapTwiceOrMore == False and overlap != -1: #cycles share exactly one node
        #create vector
        print("Create Vector, one overlap")
        vector = getCycleCycleVector(cycles[0],cycles[1],graph,overlap) #find vector using two cycle method
        print(cycles)
    elif overlapTwiceOrMore == True: #if two odd cycles are given, this should never occur!
        print("\n\nCreate Vector 2, many overlap (SHOULD NEVER OCCUR)\n\n")
        vector = "NOOOOOOOOOOOOOO"
        warning.warn("Even cycle should be present but undetected")
    else:  #if path is required to connect both cycles, work in progress
        print("Find path and create Vector 3, no overlap")
        warning.warn("Case Three Call")
        tempGraph = copy.deepcopy(graph)
        for m in range(len(cycles)):
            for n in range(len(cycles[m])):
                tempGraph[cycles[m][n],cycles[m][n-1]] = 0
                tempGraph[cycles[m][n-1],cycles[m][n]] = 0
        print(tempGraph)
        connection = findLengthConnection(cycles,tempGraph)
        print("Connection")
        print(connection)
        vector = getCyclePathCycle(cycles[0],cycles[1],connection,graph)
    return vector

#Test some more
def findLengthConnection(cycles,graph):
    path = []
    for i in range(len(cycles[0])):
        for j in range(len(cycles[1])):
            if graph[cycles[0][i],cycles[1][j]] == 1:
                path = [cycles[0][i],cycles[1][j]]
            if len(path) > 1:
                break
        if len(path) > 1:
                break
    if len(path) <= 1:
        for i in range(len(cycles)):
            for j in range(len(cycles[i])):
                path = findSingleCyclePath(graph,cycles[0][i],cycles[1][j])
                if len(path) > 1:
                    break
            if len(path) > 1:
                    break
    if len(path) <= 1:
        warning.warn("Error, graph not connected!")
    return path

def convertTetrahedraToEdges(tetrahedraList):
    print("Get list of edges")
    max = 0
    edgeList = []
    spot = 0
    for i in range(len(tetrahedraList)):
        while len(tetrahedraList[i]) > spot:
            a = tetrahedraList[i][spot]
            b = tetrahedraList[i][spot+1]
            c = tetrahedraList[i][spot+2]
            d = tetrahedraList[i][spot+3]
            spot = spot+4
            edgeList.append([a-1,b-1])
            edgeList.append([a-1,c-1])
            edgeList.append([a-1,d-1])
            edgeList.append([b-1,c-1])
            edgeList.append([b-1,d-1])
            edgeList.append([c-1,d-1])
            if d > max:
                max = d
    print("Create Graph from list")
    graph = np.zeros([max,max])
    for i in range(len(edgeList)):
       graph[edgeList[i][0],edgeList[i][1]] = 1
       graph[edgeList[i][1],edgeList[i][0]] = 1
    print(graph)
    return graph

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
    graph = convertTetrahedraToEdges(tetrahedron)
    sT = buildSpanningTree(graph)
    dendroid =  np.subtract(graph,sT)
    cycles = findOddCyclesFromSpanningTree(graph,dendroid,sT)
    cycles2 = [cycles,0]
    basis = buildBasisWithEvenCycles(graph,dendroid,sT,cycles2)
    print("Basis")
    for i in range(len(basis)):
        print(basis[i])

main()