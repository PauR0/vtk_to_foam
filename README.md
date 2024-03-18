# vtk_to_foam
A python package for conversion from vtk MultiBlock to blockMeshDict file for openFOAM. Currently, only hexahedral meshes are supported.

Following the format in which vtk/pyvista read OpenFOAM cases (Checkout the [pyvista docs](https://docs.pyvista.org/version/stable/examples/99-advanced/openfoam-example.html), those are great!), the mesh_to_foam module expects vtk MutliBlock containing two blocks, boundary and internalMesh. The block named boundary has to be itself a multiblock comprised of Polydata defining the boundary patches. For the internalMesh block only vtk UnstructuredGrid built of hexahedral elements is currently supported.
