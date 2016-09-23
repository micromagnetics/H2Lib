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

print "COORDS:", coords
print "CELLS:", cells
x = lindholm.Lindholm()
x.geometry_from_array(coords, cells)
x.setup()
