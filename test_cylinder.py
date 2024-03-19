
import sys

import pyvista as pv

from utils import change_cell_orientation
from mesh_to_foam import mesh_to_blockmeshdict


cyl = pv.CylinderStructured(radius=[1.0, .9, .8, .7], height=5, direction=[0.0,0.0,1.0]).cast_to_unstructured_grid()
cyl = change_cell_orientation(cyl) #We need to change orientation for OpenFOAM compatibility.
#cyl.plot()
cyl_boundary = cyl.extract_surface()

boundary = pv.MultiBlock()
#Let's define some patches! The last extract_surface ensures all the patches are PolyDatas
ids = cyl_boundary.points[:, 2] > 2.4 #Top lid
boundary['top lid'] = cyl_boundary.extract_points(ids, adjacent_cells=False, include_cells=True).extract_surface()

ids = cyl_boundary.points[:, 2] < -2.4 #Bottom lid
boundary['bottom lid'] = cyl_boundary.extract_points(ids, adjacent_cells=False, include_cells=True).extract_surface()

ids = (cyl_boundary.points[:, 2] > -2.4) & (cyl_boundary.points[:, 2] < 2.4) #Wall
boundary['wall'] = cyl_boundary.extract_points(ids, adjacent_cells=True, include_cells=True).extract_surface()

mesh = pv.MultiBlock()
mesh['internalMesh'] = cyl
mesh['boundary'] = boundary

if len(sys.argv) > 1:
    if sys.argv[1].lower() in ['--save', '-s']:
        mesh.save('cylinder.vtm')
else:
    #Let's visualize what we just built!
    p = pv.Plotter(shape=(1, 2))
    p.subplot(0, 0)
    p.add_mesh(mesh['internalMesh'], color='w', show_edges=True)
    p.show_bounds()
    p.subplot(0, 1)
    p.add_mesh(mesh['boundary'], multi_colors=True)
    p.add_axes()
    p.link_views()
    p.show()

bmd = mesh_to_blockmeshdict(mesh)

#print(bmd)
