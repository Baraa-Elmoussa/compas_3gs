from collections.abc import Sequence
from itertools import groupby

from compas.datastructures import Graph
from compas.datastructures import Mesh
from compas.datastructures import VolMesh
from compas.geometry import Frame
from compas.geometry import Scale
from compas.geometry import centroid_points
from compas.geometry import convex_hull
from compas_3gs.operations import cell_merge_coplanar_adjacent_faces

__all__ = [
    "point_reflection",
    "volmesh_vertex_lift",
    "volmesh_vertex_merge",
    "volmesh_halfface_pinch",
    "volmesh_halfface_extrude",
    "volmesh_halfface_subdivide",
    "volmesh_merge_adjacent_halffaces",
    "volmesh_cell_split_vertex",
    "volmesh_cell_subdivide_barycentric",
]


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   volmesh global operations
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


def point_reflection(datastructure: Graph | Mesh | VolMesh) -> None:
    """Reflect a datastructure through the centroid of its vertices, in place.

    Every point `p` moves to `2c - p`, with `c` the average of the vertex (or node) coordinates.

    Parameters
    ----------
    datastructure
        A graph, mesh or volmesh.

    Notes
    -----
    Reflecting every point keeps the direction of each face normal, so the normals stay parallel
    to the edges of the form diagram, while the cells of a volmesh (and a closed mesh) turn inside
    out: their faces point inward. This is the switch between compression and tension, see the
    attribute `convention` of a force diagram.

    References
    ----------
    - <https://en.wikipedia.org/wiki/Point_reflection>

    """
    if isinstance(datastructure, Graph):
        center = centroid_points(datastructure.nodes_attributes("xyz"))
    else:
        center = centroid_points(datastructure.vertices_attributes("xyz"))

    datastructure.transform(Scale.from_factors([-1, -1, -1], frame=Frame(center, [1, 0, 0], [0, 1, 0])))


def reverse_cycle_directions(volmesh):
    raise NotImplementedError


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   volmesh vertex operations
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


def volmesh_vertex_lift(volmesh: VolMesh, vkey: int, target_xyz: Sequence[float], hfkeys: list[int]) -> VolMesh:
    """Duplicates and lifts a vertex, then creates cells using the vertex halffaces and the duplicated vertex.

    Parameters
    ----------
    volmesh
        A volmesh object representing a polyhedral force diagram.
    vkey
        The key of the vertex to lift.
    target_xyz
        Target xyz coordinates of the lifted vertex.
    hfkeys
        List of halfface keys to create new cells from.

    Returns
    -------
    VolMesh
        The volmesh, with one new cell per halfface.

    Raises
    ------
    Exception
        If the vertex is not on the boundary.

    Notes
    -----
    The lifting vertex must be on the boundary of the volmesh.

    See Also
    --------
    `volmesh_halfface_pinch`

    """

    # check if vertex is interior ----------------------------------------------
    if not volmesh.is_vertex_on_boundary(vkey):
        raise Exception("This vertex is interior.")

    # add new, lifted vertex ---------------------------------------------------
    x, y, z = target_xyz
    w = volmesh.add_vertex(x=x, y=y, z=z)

    # add cells ----------------------------------------------------------------
    for hfkey in hfkeys:
        halffaces = [volmesh.halfface_vertices(hfkey)[::-1]]
        for u, v in volmesh.halfface_halfedges(hfkey):
            halffaces.append([u, v, w])
        volmesh.add_cell(halffaces)

    # --------------------------------------------------------------------------
    return volmesh


def volmesh_vertex_merge(volmesh: VolMesh, vkeys: list[int], target_xyz: Sequence[float]) -> VolMesh:
    """Merges specified vertices to a single vertex.

    Every cell around the vertices is rebuilt with the new vertex in place of the merged ones.
    The attributes of the cells, and of the faces and edges that remain, are kept.

    Parameters
    ----------
    volmesh
        A volmesh object representing a polyhedral force diagram.
    vkeys
        List of vertices to merge.
    target_xyz
        Target xyz coordinates of the vertices to be merged.

    Returns
    -------
    VolMesh
        The volmesh, with the merged vertex.

    Raises
    ------
    Exception
        If the vertices are not a continuous chain of halfedges in every halfface they share.

    Notes
    -----
    As implemented, the merging vertices (in any combination) must be a continuous chain of halfedges of a halfface.
    For example, vertices 1 and 3 cannot be merged.

    ```
    0 ------ 1
    |        |
    |        |
    |        |
    3 ------ 2
    ```

    """

    # inspect vertices ---------------------------------------------------------
    hf_vkeys = {}
    for vkey in vkeys:
        for hfkey in volmesh.vertex_halffaces(vkey):
            if hfkey not in hf_vkeys:
                hf_vkeys[hfkey] = []
            hf_vkeys[hfkey].append(vkey)

    for hfkey in hf_vkeys:
        all_vkeys = volmesh.halfface_vertices(hfkey)
        test_vkeys = set(hf_vkeys[hfkey])

        if len(test_vkeys) == 1:
            continue

        if len(all_vkeys) - len(test_vkeys) <= 1:
            raise Exception("Illegal merge!")

        for vkey in test_vkeys:
            anc = volmesh.halfface_vertex_ancestor(hfkey, vkey)
            des = volmesh.halfface_vertex_descendent(hfkey, vkey)
            if anc not in test_vkeys and des not in test_vkeys:
                raise Exception("Illegal merge!")

    # get vertex cells ---------------------------------------------------------
    ckeys = []
    for vkey in vkeys:
        for ckey in volmesh.vertex_cells(vkey):
            if ckey not in ckeys:
                ckeys.append(ckey)

    x, y, z = target_xyz
    new_vkey = volmesh.add_vertex(x=x, y=y, z=z)

    # construct new cell halffaces ---------------------------------------------
    cells = {}
    cell_attributes = {}
    face_attributes = {}
    edge_attributes = {}
    for ckey in ckeys:
        cell_halffaces = []
        for hfkey in volmesh.cell_faces(ckey):
            halfface = [new_vkey if key in vkeys else key for key in volmesh.halfface_vertices(hfkey)]
            halfface = [key for key, _ in groupby(halfface)]
            if len(halfface) > 1 and halfface[0] == halfface[-1]:
                halfface = halfface[:-1]
            if len(set(halfface)) > 2:
                cell_halffaces.append(halfface)
                face_attributes[frozenset(halfface)] = volmesh.face_attributes(hfkey)
        cells[ckey] = cell_halffaces
        cell_attributes[ckey] = volmesh.cell_attributes(ckey)
        for edge in volmesh.cell_edges(ckey):
            u, v = [new_vkey if key in vkeys else key for key in edge]
            if u != v:
                edge_attributes[frozenset((u, v))] = volmesh.edge_attributes(edge)

    for ckey in ckeys:
        volmesh.delete_cell(ckey)

    # add cells and copy the attributes back -----------------------------------
    for ckey in cells:
        volmesh.add_cell(cells[ckey], ckey=ckey)
        volmesh.cell_attributes(ckey, names=cell_attributes[ckey].keys(), values=cell_attributes[ckey].values())
        for hfkey in volmesh.cell_faces(ckey):
            attr = face_attributes.get(frozenset(volmesh.halfface_vertices(hfkey)), {})
            for name, value in attr.items():
                volmesh.face_attribute(hfkey, name, value)
        for edge in volmesh.cell_edges(ckey):
            attr = edge_attributes.get(frozenset(edge), {})
            for name, value in attr.items():
                volmesh.edge_attribute(edge, name, value)

    volmesh.remove_unused_vertices()

    # --------------------------------------------------------------------------
    return volmesh


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   volmesh edge operations
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


# def _volmesh_cull_duplicate_edges(volmesh):
#     new_edges = {vkey: {} for vkey in volmesh.vertex}
#     seen      = set()
#     for u, v in volmesh.edges():
#         edge = frozenset([u, v])
#         if edge not in seen:
#             new_edges[u][v] = {}
#             seen.add(edge)
#     volmesh.edge = new_edges


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   volmesh halfface operations
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


def volmesh_halfface_extrude(volmesh, hfkey):
    raise NotImplementedError


def volmesh_halfface_subdivide(volmesh, hfkey):
    raise NotImplementedError


def volmesh_halfface_pinch(volmesh: VolMesh, hfkey: int, xyz: Sequence[float]) -> VolMesh:
    """Add a pyramid cell on a boundary halfface, with its apex at a new vertex.

    Parameters
    ----------
    volmesh
        A volmesh object representing a polyhedral force diagram.
    hfkey
        The boundary halfface to pinch.
    xyz
        XYZ coordinates of the apex.

    Returns
    -------
    VolMesh
        The volmesh, with the new cell.

    """

    x, y, z = xyz
    w = volmesh.add_vertex(x=x, y=y, z=z)
    hf_vkeys = volmesh.halfface_vertices(hfkey)
    cell_halffaces = [hf_vkeys[::-1]]
    halfedges = volmesh.halfface_halfedges(hfkey)
    for u, v in halfedges:
        cell_halffaces.append([u, v, w])
    volmesh.add_cell(cell_halffaces)
    return volmesh


def volmesh_merge_adjacent_halffaces(volmesh: VolMesh, hfkeys: list[int]) -> VolMesh:
    """Add a cell on a chain of adjacent boundary halffaces: the convex hull of their vertices.

    Coplanar faces of the hull are merged. The new cell shares the given halffaces with the volmesh.

    Parameters
    ----------
    volmesh
        A volmesh object representing a polyhedral force diagram.
    hfkeys
        The boundary halffaces, each one adjacent to at least one other.

    Returns
    -------
    VolMesh
        The volmesh, with the new cell.

    Raises
    ------
    ValueError
        If a halfface is interior, or if the halffaces are not chained.

    """

    # check halffaces ----------------------------------------------------------
    for hfkey in hfkeys:
        if not volmesh.is_halfface_on_boundary(hfkey):
            raise ValueError("Halfface {} is interior.".format(hfkey))
    if not _are_halffaces_chained(volmesh, hfkeys):
        raise ValueError("These halffaces are not chained.")
    # --------------------------------------------------------------------------
    halffaces = [volmesh.halfface_vertices(hfkey) for hfkey in hfkeys]
    # --------------------------------------------------------------------------
    vkeys = set()
    for hfkey in hfkeys:
        for key in volmesh.halfface_vertices(hfkey):
            vkeys.add(key)

    vkeys = list(vkeys)
    points = [volmesh.vertex_coordinates(vkey) for vkey in vkeys]

    faces_by_index = convex_hull(points)
    faces_by_vkeys = []
    for face in faces_by_index:
        faces_by_vkeys.append([vkeys[index] for index in face])

    # make temp cell mesh ------------------------------------------------------
    cell = Mesh()
    for i in range(len(vkeys)):
        key = vkeys[i]
        x, y, z = points[i]
        cell.add_vertex(key=key, x=x, y=y, z=z)
    for face in faces_by_vkeys:
        cell.add_face(face)

    # merge coplanar faces -----------------------------------------------------
    cell_merge_coplanar_adjacent_faces(cell)

    # get correct direction of faces -------------------------------------------
    # a face of the hull may list the vertices of halffaces[0] from another start vertex
    faces = [cell.face_vertices(fkey) for fkey in cell.faces()]
    first = halffaces[0]
    rotations = [first[i:] + first[:i] for i in range(len(first))]
    if any(face in rotations for face in faces):
        new_faces = [face[::-1] for face in faces]
    else:
        new_faces = faces

    volmesh.add_cell(new_faces)

    return volmesh


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   volmesh cell operations
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


def volmesh_cell_split_vertex(volmesh: VolMesh, vkey: int, hfkey1: int, hfkey2: int) -> None:
    """Split a vertex of a single-cell volmesh between two of its halffaces.

    Raises
    ------
    NotImplementedError
        Always: the method relied on an EGI stored per cell in `c_data` (compas 0.x).

    """
    raise NotImplementedError("volmesh_cell_split_vertex relied on the cell EGI stored in `c_data` (compas 0.x) and is not available.")


def volmesh_cell_subdivide_barycentric(volmesh: VolMesh, ckey: int) -> list[int]:
    """Subdivide a cell of a volmesh into pyramids, one per face, with their apex at the cell center.

    The attributes of the faces and edges of the cell are kept; the new cells get the default
    cell attributes.

    Parameters
    ----------
    volmesh
        A volmesh object representing a polyhedral force diagram.
    ckey
        The identifier of the cell.

    Returns
    -------
    list[int]
        The identifiers of the new cells.

    """
    new_ckeys = []
    x, y, z = volmesh.cell_center(ckey)

    halffaces = [volmesh.halfface_vertices(hfkey) for hfkey in volmesh.cell_faces(ckey)]
    face_attributes = {frozenset(volmesh.halfface_vertices(hfkey)): volmesh.face_attributes(hfkey) for hfkey in volmesh.cell_faces(ckey)}
    edge_attributes = {frozenset(edge): volmesh.edge_attributes(edge) for edge in volmesh.cell_edges(ckey)}

    volmesh.delete_cell(ckey)

    w = volmesh.add_vertex(x=x, y=y, z=z)

    for halfface in halffaces:
        cell_halffaces = [halfface]

        halfedges = zip(halfface, halfface[1:] + halfface[:1])

        for u, v in halfedges:
            cell_halffaces.append([w, v, u])

        new_ckeys.append(volmesh.add_cell(cell_halffaces))

    # copy the attributes back -------------------------------------------------
    for new_ckey in new_ckeys:
        for hfkey in volmesh.cell_faces(new_ckey):
            attr = face_attributes.get(frozenset(volmesh.halfface_vertices(hfkey)), {})
            for name, value in attr.items():
                volmesh.face_attribute(hfkey, name, value)
        for edge in volmesh.cell_edges(new_ckey):
            attr = edge_attributes.get(frozenset(edge), {})
            for name, value in attr.items():
                volmesh.edge_attribute(edge, name, value)

    return new_ckeys


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   volmesh operation helpers
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


def _are_halffaces_chained(volmesh: VolMesh, hfkeys: list[int]) -> bool:
    chained = set()

    for hfkey_1 in hfkeys:
        nbrs = []

        # edge_halffaces lists the halffaces with the halfedge in the given direction only
        for u, v in volmesh.halfface_halfedges(hfkey_1):
            nbrs += volmesh.edge_halffaces((u, v)) + volmesh.edge_halffaces((v, u))

        for hfkey_2 in hfkeys:
            if hfkey_2 != hfkey_1 and hfkey_2 in nbrs:
                chained.add(hfkey_1)
                chained.add(hfkey_2)

    if len(hfkeys) == 1 or len(chained) == len(hfkeys):
        return True

    return False


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
