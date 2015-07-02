__author__ = 'Owner'

def checkTriangles(edge1,edge2,edge3,edge4,edge5,edge6):
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


def checkLegalTetrahedron(edge1,edge2,edge3,edge4,edge5,edge6):
        legalTetrahedron=True
        det = (-2)*(edge1**4)*(edge6**2)-(2)*(edge1**2)*(edge2**2)*(edge4**2)+(2)*(edge1**2)*(edge2**2)*(edge5**2)+(2)*(edge1**2)*(edge2**2)*(edge6**2)+(2)*(edge1**2)*(edge3**2)*(edge4**2)-(2)*(edge1**2)*(edge3**2)*(edge5**2)+(2)*(edge1**2)*(edge3**2)*(edge6**2)+(2)*(edge1**2)*(edge4**2)*(edge6**2)+(2)*(edge1**2)*(edge5**2)*(edge6**2)-(2)*(edge1**2)*(edge6**4)-(2)*(edge2**4)*(edge5**2)+(2)*(edge2**2)*(edge3**2)*(edge4**2)+(2)*(edge2**2)*(edge3**2)*(edge5**2)-(2)*(edge2**2)*(edge3**2)*(edge6**2)+(2)*(edge2**2)*(edge4**2)*(edge5**2)-(2)*(edge2**2)*(edge5**4)+(2)*(edge2**2)*(edge5**2)*(edge6**2)-(2)*(edge3**4)*(edge4**2)-(2)*(edge3**2)*(edge4**4)+(2)*(edge3**2)*(edge4**2)*(edge5**2)+(2)*(edge3**2)*(edge4**2)*(edge6**2)-(2)*(edge4**2)*(edge5**2)*(edge6**2)
        if det<=0:
            legalTetrahedron=False
        return legalTetrahedron

edge12=3.767999428
edge13=8.145635447
edge14=4.017298606
edge15=6.38192561
edge23=5.752376321
edge24=4.281292903
edge25=4.439533053
edge34=4.932513323
edge35=3.346511184
edge45=4.892056976

print('legal Triangles')
print(checkTriangles(edge12,edge13,edge14,edge23,edge24,edge34))
print(checkTriangles(edge12,edge13,edge15,edge23,edge25,edge35))
print(checkTriangles(edge12,edge14,edge15,edge24,edge25,edge45))
print(checkTriangles(edge13,edge14,edge15,edge34,edge35,edge45))
print(checkTriangles(edge23,edge24,edge25,edge34,edge35,edge45))

print('legal tetrahedron')
print(checkLegalTetrahedron(edge12,edge13,edge14,edge23,edge24,edge34))
print(checkLegalTetrahedron(edge12,edge13,edge15,edge23,edge25,edge35))
print(checkLegalTetrahedron(edge12,edge14,edge15,edge24,edge25,edge45))
print(checkLegalTetrahedron(edge13,edge14,edge15,edge34,edge35,edge45))
print(checkLegalTetrahedron(edge23,edge24,edge25,edge34,edge35,edge45))