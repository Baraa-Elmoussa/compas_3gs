from collections.abc import Sequence

from compas.geometry import distance_point_point
from compas.geometry import dot_vectors
from compas.geometry import intersection_line_plane
from compas.geometry import midpoint_point_point
from compas.tolerance import TOL
from compas_3gs.diagrams import EGI
from compas_3gs.diagrams import Cell

__all__ = ["cell_split_indet_face_vertices", "cell_collapse_short_edge", "cell_relocate_face", "cell_face_subdivide_barycentric", "cell_merge_coplanar_adjacent_faces"]


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   cell vertices
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


def cell_split_indet_face_vertices(cell: Cell, fkey: int) -> Cell:
    """Split all indeterminate vertices of a cell face.

    Parameters
    ----------
    cell
        Cell as a mesh object.
    fkey
        Identifier of the face.

    Returns
    -------
    Cell
        Updated cell.

    Notes
    -----
    An indeterminate vertex is defined as a vertex with degree or valency of 3 or greater.

    """
    egi = cell.dual(cls=EGI)

    f_vkeys = cell.face_vertices(fkey)
    egi_nbr_vkeys = egi.vertex_neighbors(fkey)

    for egi_fkey in f_vkeys:
        fkeys = egi.face_vertices(egi_fkey)
        i = fkeys.index(fkey)
        fkeys = fkeys[i:] + fkeys[:i]

        egi_face_vertices = [key for key in fkeys if key not in egi_nbr_vkeys + [fkey]]

        vkey_del = egi_fkey
        x, y, z = cell.vertex_coordinates(vkey_del)

        for vkey in egi_face_vertices:
            f, g = egi.split_face(vkey_del, fkey, vkey)

            cell.delete_vertex(vkey_del)

            cell.add_vertex(key=f, x=x, y=y, z=z)
            cell.add_vertex(key=g, x=x, y=y, z=z)

            for new_fkey in fkeys:
                new_vkeys = egi.vertex_faces(new_fkey, ordered=True)
                cell.add_face(new_vkeys, fkey=new_fkey)
            vkey_del = g

    return cell


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   cell edges
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


def cell_collapse_short_edge(cell: Cell, u: int, v: int, min_length: float = 0.1) -> Cell:
    """Collapse a short edge of a cell geometrically.

    If the edge is shorter than `min_length`, both vertices move to its midpoint,
    respecting their fixity. The topology of the cell is unchanged.

    Parameters
    ----------
    cell
        Cell as a mesh object.
    u
        The key of the start vertex.
    v
        The key of the end vertex.
    min_length
        Minimum length of edges to be collapsed.

    Returns
    -------
    Cell
        Updated cell.

    """
    sp = cell.vertex_coordinates(u)
    ep = cell.vertex_coordinates(v)
    dist = distance_point_point(sp, ep)

    if dist < min_length:
        mp = midpoint_point_point(sp, ep)
        cell.vertex_update_xyz(u, mp)
        cell.vertex_update_xyz(v, mp)

    return cell


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   cell face operations
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


def cell_cull_zero_faces(cell: Cell) -> None:
    pass


def cell_relocate_face(cell: Cell, fkey: int, xyz: Sequence[float], normal: Sequence[float]) -> Cell:
    """Relocate the face of a cell to a target plane.

    Parameters
    ----------
    cell
        Cell as a mesh object.
    fkey
        Identifier of the face.
    xyz
        XYZ coordinates of a point on the target plane.
    normal
        Normal vector of the target plane.

    Returns
    -------
    Cell
        Updated cell.

    """

    cell_split_indet_face_vertices(cell, fkey)

    vkeys = cell.face_vertices(fkey)

    # new target plane for the face
    plane = (xyz, normal)

    # neighboring edges
    edges = {}
    for u in vkeys:
        for v in cell.vertex_neighbors(u):
            if v not in vkeys:
                edges[u] = v

    for u in edges:
        line = cell.edge_coordinates((u, edges[u]))
        it = intersection_line_plane(line, plane)
        cell.vertex_update_xyz(u, it, constrained=False)

    return cell


def cell_merge_coplanar_adjacent_faces(cell: Cell, tol: float | None = None) -> Cell:
    """Merge adjacent faces of a cell that are coplanar.

    Each face is merged with the neighbors for which `1 - dot(n1, n2) < tol`, with `n1` and `n2`
    the unit normals of the two faces. A merged face gets a new key.

    Parameters
    ----------
    cell
        Cell as a mesh object.
    tol
        Tolerance on `1 - dot(n1, n2)`. Default is `TOL.approximation`.

    Returns
    -------
    Cell
        Updated cell.

    """
    tol = TOL.approximation if tol is None else tol

    initial_faces = list(cell.faces())
    current_faces = list(cell.faces())

    for fkey in initial_faces:
        if fkey in current_faces:
            normal = cell.face_normal(fkey)

            # a neighbor that shares several edges with the face is listed once per edge
            nbr_fkeys = list(dict.fromkeys(cell.face_neighbors(fkey)))

            faces_to_merge = []
            for nbr_fkey in nbr_fkeys:
                nbr_normal = cell.face_normal(nbr_fkey)
                dot = dot_vectors(normal, nbr_normal)
                if 1 - dot < tol:
                    faces_to_merge.append(nbr_fkey)

            if faces_to_merge:
                current_faces.remove(fkey)
                key = fkey
                # merge two faces at a time
                for nbr_fkey in faces_to_merge:
                    key = cell.merge_faces([key, nbr_fkey])
                    current_faces.remove(nbr_fkey)
                current_faces.append(key)

    return cell


def cell_face_subdivide_barycentric(cell: Cell, fkey: int) -> Cell:
    """Subdivide a face of a cell into a fan of triangles around its centroid.

    Parameters
    ----------
    cell
        Cell as a mesh object.
    fkey
        Identifier of the face.

    Returns
    -------
    Cell
        Updated cell.

    """
    cell.insert_vertex(fkey, xyz=cell.face_centroid(fkey))
    return cell


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   Main
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


if __name__ == "__main__":
    pass
