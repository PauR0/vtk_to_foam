# vtk_to_foam
A python package for conversion from vtk MultiBlock to blockMeshDict file for openFOAM. Currently, **only hexahedral** meshes are supported.

Following the format in which vtk/pyvista read OpenFOAM cases (Checkout the [pyvista docs](https://docs.pyvista.org/version/stable/examples/99-advanced/openfoam-example.html), those are great!), the mesh_to_foam module expects vtk MutliBlock containing two blocks, boundary and internalMesh. The block named boundary has to be itself a multiblock comprised of Polydata defining the boundary patches. For the internalMesh block only vtk UnstructuredGrid built of hexahedral elements is currently supported.

Tested using OpenFOAM v10 and python3.8 (with pyvista==0.43.4)

Assuming the availability of a testing OpenFOAM case directory at OF_case, the following commands should work.

    python test_cylinder.py -s
    ./mesh_to_foam.py cylinder.vtm -c OF_case/system/
    cd OF_case && blockMesh
