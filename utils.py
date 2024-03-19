
import os
import sys
import re
from uuid import uuid4
from shutil import copytree


import numpy as np
import pyvista as pv
from scipy.spatial import KDTree
from scipy.spatial.transform import Rotation


path = os.path.realpath(__file__)
template_dir = os.path.join(os.path.dirname(path), 'template_case')


def normalize(v, o=2):
    """
    Normalize an array. Using numpy linalg's norm.
    If the norm is zero an array or zeros is returned.

    Arguments:
    -----------

        v : np.ndarray
            The array to be normalized.

        o : any
            See numpy linalg norms accepted values for argument ord.

    Returns:
    ---------
        : np.ndarray
            The normalized array.
    """
    n = np.linalg.norm(v, ord=o)
    if ord == 0.0:
        return v*0
    return v/n
#

def get_uid(n=10):
    """ Get an unique hex ID of length n
    """
    return uuid4().hex[:n]
#

def find_and_rep(text, rep):
    """
    Find and replace substrings using regex. The replacements are
    described by means of the rep variable, a dictionary where
    the key is the string to be found and the value is the string
    to be placed instead.



    Arguments:
    -----------
        text : str
            The original text that is to be modified.

        rep : dict
            The strings to find as keys and the ones to replace
            as values.


    Returns:
    ---------
        text : str
            The modified string


    Example:
    ---------
    text = "Crucifixion? No, freedom for me"
    rep = {
        'No' : 'Yes',
        'freedom' : 'crucifixion'
    }
    new_text = regex_find_and_rep(text, rep)
    new_text
    'Crucifixion? Yes, crucifixion for me'

    """

    rep = dict((re.escape(k), v) for k, v in rep.items())
    pattern = re.compile("|".join(rep.keys()))
    text = pattern.sub(lambda m: rep[re.escape(m.group(0))], text)

    return text
#

def adapt_file(source, rep, target=None, logg=None):
    """
    This function modifies a base source file by means of find_and_rep function,
    and, if target is provided saves the modified to disk at saved path. The function
    returns the a string containing the modified file

    WARNNG: This function does not perform any check on whether the target file
    exists or not and may overwrite some useful data.


    Arguments:
    -----------

        source : str, optional
            The base file.

        rep : dict
            The dictionary containing the find and replace str.
            If the dictionary is empty, nothing is modified.

        target : str, opt
            Default None. The file name of the target file.

        logg : logger
            Default None. The logger to logg the result, if none is passed
            it uses the standard output.

    Returns:
    ---------

        out_str

    """

    with open(source, 'r') as f_in:
        in_str = f_in.read()

    if rep:
        out_str = find_and_rep(text=in_str, rep=rep)
    else:
        out_str = in_str

    if target is not None:
        with open(target, 'w') as f_out:
            f_out.write(out_str)
        if logg:
            logg.info(f"  Created .........................{target}")
        else:
            print(f"  Created .........................{target}")

    return out_str
#

def get_point_in_mesh(mesh):
    """
    Get a point located inside the simulation mesh. The function randomly selects a point
    on the mesh and traces the ray using its normal. Then the considered point is the center
    between the initial and the intersected by the ray.

    Arguments:
    ------------

        openfoam_case_path : str
            The openfoam directory containing the file constant/triSurface/aorta_mesh.stl

    """

    voxels = pv.voxelize(mesh, density=0.1, check_surface=True)
    vox_surf = voxels.extract_surface(pass_pointid=True)

    for i in range(voxels.n_points):
        if i not in vox_surf['vtkOriginalCellIds']:
            return voxels.points[i]

    print("WARNING: could not find point inside the mesh. Returning None.")
    return None
#

def make_0_dir_backup(openfoam_case_path, logg=None):
    """
    Make a backup 0.org directory with the initial boundary conditions.

    Arguments:
    -----------

        openfoam_case_path : string
                Destionation path where the OpenFOAM case is

        logg : logger
            To logg.

    """
    source = f"{openfoam_case_path}/0"
    dest = f"{openfoam_case_path}/0.orig"

    copytree(source, dest)

    if logg:
        logg.info(f"  Created .........................{dest}")
    else:
        print(f"  Created .........................{dest}")
#

def touch(fname):
    """

    Function equivalent to unix touch

        Arguments:
    -----------

        fname : string
                The file to touch

    """
    try:
        os.utime(fname, None)
    except OSError:
        open(fname, 'a').close()
#

def make_dot_foam_file(openfoam_case_path):
    """

    Function to create the case.foam file inside the
    case.

    Arguments:
    -----------

        openfoam_case_path : string
                Destionation path where the OpenFOAM case is

    """

    of_file = f"{openfoam_case_path}/case.foam"
    touch(of_file)

    return of_file
#

def points_to_of_str(pts):
    """
    Convert a numpy array with dimensions (n, 3) to openFOAM string format
    which is:

    (
    (x1, y1, z1)
    ...
    (xn, yn, zn)
    );

    Arguments:
    -------------

        pts : np.ndarray, (n, 3)
            The  array of points

    Returns:
    ----------

        pts_str : string
            The string of points formatted in the OpenFOAM syntax.
    """

    pts_str = np.array2string(pts, separator=' ', threshold=sys.maxsize)[1:-1].replace('[','(').replace(']',')')
    return f"(\n {pts_str}\n);\n"
#

def make_vertices_for_blockmesh(mesh, n=None, ext=1.1):
    """
    Return the vertices list for the blockMeshDict from either a
    pv.PolyData or pv.MultiBlock. If a PolyData has been passed,
    the vertices are those of it bounding box. If a multiblock has been
    passed, it is assumed that it contains the internalMesh key, corresponding
    to an unstructured grid built of hexahedral cells.

    Arguments:
    -------------

        mesh : pv.PolyData | pyvista.MultiBlock
            The mesh to be simulated either volumetric
            or the surface.

        n : np.array (3,), opt.
            Default None. If the n vector is passed, the blockMesh is rotated
            to align the Z axis with it.

        ext : float, optional.
            An extension factor used to slightly increase the bounding box
            to prevent element intersection with the boundary patches. Defaulting
            to 1.1.


    Return:
    -------

        vertices : np.ndarray (n, 3)
            The vertices array.
    """

    verts = None

    if isinstance(mesh, pv.MultiBlock):
        verts = mesh['internalMesh'].points

    elif isinstance(mesh, pv.PolyData):

        #Scale to prevent intersection between elements and patches
        points = (mesh.points - mesh.center) * ext + mesh.center

        if n is not None:
            e3 = np.array([0,0,1])
            rv = normalize(np.cross(mesh.inlet_normal, e3)) * np.arccos(e3.dot(mesh.inlet_normal))
            rot = Rotation.from_rotvec(rv)
            pts_tr = rot.apply(points)
            (x_min, y_min, z_min), (x_max, y_max, z_max) = pts_tr.min(axis=0), pts_tr.max(axis=0)

        else:
            (x_min, y_min, z_min), (x_max, y_max, z_max) = points.min(axis=0), points.max(axis=0)

        verts = np.array(((x_min, y_min, z_min),
                            (x_max, y_min, z_min),
                            (x_max, y_max, z_min),
                            (x_min, y_max, z_min),
                            (x_min, y_min, z_max),
                            (x_max, y_min, z_max),
                            (x_max, y_max, z_max),
                            (x_min, y_max, z_max)))

        if n is not None:
            verts = rot.apply(verts, inverse=True)

    else:
        print(f"ERROR: make_vertices_for_blockmesh only supports pv.PolyData and pv.MultiBlock, and type {mesh.__class__} was passed.")

    return verts
#

def make_blocks_for_blockmesh(mesh, n):
    """
    A function to build hex blocks for block mesh.

    Arguments:
    -----------

        mesh : AortaMesh | pv.MultiBlock
            The aorta mesh

        n : tuple(int)
            Int or tuple of ints (n1, n2, n3) with the number
            of division for each block side.

    Returns:
    ---------

        blocks : str
            The string of blocks to be written on the blockMeshDict.

    """

    if isinstance(n, int):
        n = [n]*3

    if isinstance(mesh, pv.PolyData):
        return f"(\nhex (0 1 2 3 4 5 6 7) ({n[0]} {n[1]} {n[2]}) simpleGrading (1 1 1)\n);\n"

    if isinstance(mesh, pv.MultiBlock):
        cell_types = np.unique( mesh['internalMesh'].celltypes)
        if not np.in1d(cell_types, [pv.CellType.WEDGE, pv.CellType.HEXAHEDRON]).all():
            print(f"ERROR: Some cells are not hexahedra nor wedge, this is not currently supported. Celltypes found: {cell_types}")
            return

        hexs  = []
        if pv.CellType.WEDGE in cell_types:
            hexs += [f"hex ({c[0]} {c[1]} {c[2]} {c[0]} {c[3]} {c[4]} {c[5]} {c[3]}) (1 1 1) simpleGrading (1 1 1)"
                                                        for c in mesh['internalMesh'].cells_dict[pv.CellType.WEDGE]]
        if pv.CellType.HEXAHEDRON in cell_types:
            hexs += [f"hex ({str(c)[1:-1]}) (1 1 1) simpleGrading (1 1 1)"
                                                        for c in mesh['internalMesh'].cells_dict[pv.CellType.HEXAHEDRON]]

        hex_str = "\n".join(hexs)
        return f"(\n{hex_str}\n);\n"
    else:
        print(f"ERROR: make_blocks_for_blockmesh only supports pv.PolyData and pv.MultiBlock data, and type {mesh.__class__} was passed.")
        return
#

def extract_faces(patch, points, kdt=None):
    """
    This function extract the face lists for a given patch.

    TODO: Check to speed this up using "vtkOriginalPointIds".

    Arguments:
    ------------
        patches : pv.PolyData
            A given boundary patch.

        points : np.ndarray
            The points of the whole mesh, not only the patch.

    Returns:
    -----------
        faces : list[int]
            The faces list with the point ids building each face.
    """

    if kdt is None:
        kdt = KDTree(points)

    i, n = 0, 0
    faces = []
    while i < patch.n_faces_strict:
        n_vert = patch.faces[n]
        f = patch.faces[n+1:n+n_vert+1]
        face = kdt.query(patch.points[f])[1]
        if len(face) == 3: #CAUTION when using non thetrahedral meshes, this may fail....
            face = np.concatenate([face, [face[0]]])
        faces.append(face)
        #Update
        n += n_vert+1
        i += 1

    return faces
#

def make_patches_for_block_mesh(mesh):
    """
    A function to build face patches for blockMesh.

    Arguments:
    -----------

        mesh : pv.MultiBlock
            The multiblock mesh with 'boundary' and 'internalMesh'.

        n : tuple(int)
            Int or tuple of ints (n1, n2, n3) with the number
            of division for each block side.

    Returns:
    ---------
        bounds_str : str
            The string of patches to be written on the blockMeshDict.
    """

    boundary, internalMesh = mesh['boundary'], mesh['internalMesh']

    kdt_locator = KDTree(internalMesh.points)

    bounds = []
    for bound_id, polybound in zip(boundary.keys(), boundary):
        if bound_id == 'wall':
            ptype = 'wall'
        else:
            ptype = 'patch'

        faces = extract_faces(polybound, internalMesh.points, kdt_locator)
        facestr = "\n\t\t".join([f"({str(f)[1:-1]})" for f in faces])
        bound_str = f"\t{bound_id}\n\t{{\n\t\ttype {ptype};\n\t\tfaces\n\t\t(\n\t\t{facestr}\n\t\t);\n\t}}\n"
        bounds.append(bound_str)

    bounds_str = "\n".join(bounds)
    return f"(\n{bounds_str});"
#

def handle_path(fname, save_dir=None, w=False):
    """
    Function to perform the path arrangmente of filename and
    check if file exists.

    If save_dir is None, it returns True. If save_dir is not None,
    then it concatenates

    Arguments:
    -----------

        fname : str
            The name of the file

        save_dir : str
            Path to the directory where file should be saved.

        w : bool, opt
            Default False. Whether to overwrite existing files.

    Returns:
    ------------
        fout : str
            The complete file name or an empty string if file already exists
            and overwritting is set to False.
    """


    if save_dir is not None:
        fname = os.path.join(save_dir, fname)

    if os.path.exists(fname) and not w:
        print(f"WARNING: {fname} already exists and overwrite is set to False. Nothing will be saved.")
        return ''

    return fname
#