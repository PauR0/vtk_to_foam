#! /usr/bin/env python3
import os
import argparse

import pyvista as pv

from utils import (template_dir,
                   make_vertices_for_blockmesh,
                   make_blocks_for_blockmesh,
                   make_patches_for_block_mesh,
                   points_to_of_str,
                   adapt_file)


def make_block_mesh_dict(mesh, fname=None, n=None, logg=None):
    """
    Set the parameters of the snappyHexMesh box
    to make them be equal to the aorta bounding box

    Arguments:
    ----------

        mesh : AortaMesh
            A mesh containing an aorta geometry.

        fname : str, opt.
            A filename to save the blockMeshDict.

        n : np.array (3,), opt.
            Default None. If the n vector is passed, the blockMesh is rotated
            to align the Z axis with it.

    """

    verts = make_vertices_for_blockmesh(mesh, n=n, ext=1.1)
    n = 1
    rep = {
        "**vertices**" : points_to_of_str(verts),
        "**blocks**"   : make_blocks_for_blockmesh(mesh, n),
        "**boundary**" : make_patches_for_block_mesh(mesh)
    }

    source = os.path.join(template_dir, 'system', 'blockMeshDict')
    return adapt_file(source=source, rep=rep, target=fname, logg=logg)
#


def mesh_to_blockmeshdict(mesh, save_dir=None, w=False):
    """
    This function turns a pyvista (vtk) Multiblock built of hexahedra
    and converts it to the blockMeshDict format.

    Arguments:
    ------------

        mesh : pv.MultiBlock
            The multiblock containing the 'bundary' and 'internalMesh' entries
            with the boundary patches and the hexahedral mesh.

        save_dir : str
            Path to the system directory of the OpenFOAM simulation case.

        w : bool, opt
            Default False. Whether to overwrite existing files.

    """

    fname=None
    if save_dir is not None:
        fname = os.path.join(save_dir, 'blockMeshDict')
        if os.path.exists(fname) and not w:
            print(f"WARNING: {fname} already exists and overwritting is set to False. Nothing will be written...")
            return

    return make_block_mesh_dict(mesh=mesh, fname=fname)
#

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""blockMeshDict generator for vtk MultiBlocks""")

    parser.add_argument('-w',
                        action='store_true',
                        help="""Force file overwrite if already exists.""")

    parser.add_argument('-c',
                        '--case',
                        dest='case',
                        type=str,
                        default=None,
                        help = """Path to the openfoam case. If nothing is given,
                        the current directory is used.""")

    parser.add_argument('mesh',
                    nargs=1,
                    type=str,
                    help = """The path to the vtk MultiBlock mesh.""")


    args = parser.parse_args()

    mesh = pv.read(args.mesh)

    mesh_to_blockmeshdict(mesh=mesh, save_dir=args.case, w=args.w)
#
