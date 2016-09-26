import lindholm
from dolfin import *

mesh = UnitCubeMesh(10,10,10)
bmesh = BoundaryMesh(mesh, "exterior", False)
coords = bmesh.coordinates()
cells = bmesh.cells()

File("mesh.pvd") << mesh
File("bmesh.pvd") << bmesh
#x = lindholm.Lindholm()
#x.geometry_from_file("model.msh")
#x.setup()

x = lindholm.Lindholm()
x.geometry_from_array(coords, cells)
x.setup()
print "TEST"

V = FunctionSpace(mesh, "CG", 1)
u1 = interpolate(Constant(1.0), V)

print x.matvec(u1.vector().array())
