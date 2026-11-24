
from collections.abc import Sequence
from math import atan2

from compas.datastructures import Graph

from compas.datastructures.graph.duality import graph_node_find_first_neighbor

from compas.geometry import Frame
from compas.geometry import Plane
from compas.geometry import add_vectors
from compas.geometry import angle_vectors
from compas.geometry import distance_point_point as distance
from compas.geometry import dot_vectors
from compas.geometry import cross_vectors
from compas.geometry import length_vector
from compas.geometry import normalize_vector
from compas.geometry import scale_vector
from compas.geometry import subtract_vectors

from compas.tolerance import TOL

from compas_3gs.diagrams import EGI
from compas_3gs.diagrams import Cell


__all__ = ['egi_from_vectors',
           'cell_from_egi']


def egi_from_vectors(vectordict: dict[int, Sequence[float]], origin: Sequence[float], tol: float | None = None) -> EGI:
    """Construct an egi from a set of vectors.

    Parameters
    ----------
    vectordict
        A dictionary of key-vector pairs.
    origin
        The coordinates of the centroid.
    tol
        Tolerance for detecting coincident vectors, on `1 - dot(n1, n2)` of the unit vectors.
        Default is `TOL.approximation`.

    Returns
    -------
    EGI
        A mesh object representing the egi.

    Raises
    ------
    Exception
        If there are less than four vectors, if two vectors coincide,
        or if the faces of the egi are not consistently oriented.

    Notes
    -----
    The adjacency arcs are great-circle arcs between the unit vectors, and the cross-adjacencies
    are the intersections of these arcs, both computed with vector algebra on the unit sphere.
    Two points on the unit sphere are the same if they are closer than `TOL.angular`.

    Warning
    -------
    - This algorithm does not address scenarios where multiple parallel (collinear) vectors are present.

    References
    ----------
    - Horn, B.K.P. (1984). *Extended Gaussian images*.
    - Moni, S. (1990, June). *A closed-form solution for the reconstruction
    of a convex polyhedron from its extended gaussian image.*
    - Lee, J., T. Van Mele, and P. Block (2018). *Disjointed force polyhedra.*

    """
    tol = TOL.approximation if tol is None else tol

    if len(vectordict) < 4:
        raise Exception('Four or more vectors are needed for the construction of egi.')

    egi = Graph()

    # --------------------------------------------------------------------------
    #   1. add vertices from vectors
    # --------------------------------------------------------------------------
    for vkey in vectordict:
        normal = normalize_vector(vectordict[vkey])
        x, y, z = add_vectors(normal, origin)
        egi.add_node(vkey, x=x, y=y, z=z, type='face', normal=normal, nbrs=[])

    # --------------------------------------------------------------------------
    #   2.  Identify main adjacencies
    # --------------------------------------------------------------------------
    vkey_pairs = set()

    for vkey in egi.nodes():
        v_crs_list = []

        for nbr_vkey in egi.nodes():

            if nbr_vkey is not vkey:

                n1 = egi.node_attribute(vkey, 'normal')
                n2 = egi.node_attribute(nbr_vkey, 'normal')

                # This checks if the normals are opposite ----------------------
                dot = dot_vectors(n1, n2)

                if dot > 1 - tol:
                    raise Exception("Coincident vectors detected.")

                elif dot > -1 + tol:

                    unit_crs = normalize_vector(cross_vectors(n1, n2))

                    # Check to see if any other normals are coplanar
                    index = None
                    for i, (crs, _) in enumerate(v_crs_list):
                        if distance(crs, unit_crs) < TOL.angular:
                            index = i
                            break

                    if index is None:
                        v_crs_list.append((unit_crs, nbr_vkey))

                    # If multiple arcs are coplanar, choose the closer one
                    else:
                        this_dist = distance(egi.node_coordinates(vkey),
                                             egi.node_coordinates(nbr_vkey))
                        test_dist = distance(egi.node_coordinates(vkey),
                                             egi.node_coordinates(v_crs_list[index][1]))
                        if this_dist < test_dist:
                            v_crs_list[index] = (v_crs_list[index][0], nbr_vkey)

        # Add to overall connectivity dict -------------------------------------
        for _, nbr_vkey in v_crs_list:
            pair = frozenset([vkey, nbr_vkey])
            vkey_pairs.add(pair)

    # --------------------------------------------------------------------------
    #   3.  Main adjacency arcs
    # --------------------------------------------------------------------------
    arcs = {}

    for pair in vkey_pairs:
        u, v = list(pair)
        arc = _draw_arc(egi.node_attribute(u, 'normal'),
                        egi.node_attribute(v, 'normal'),
                        origin)

        if len(arcs) == 0:
            arc_key = 0
        else:
            arc_key = max(arcs) + 1
        arcs[arc_key] = {'arc': arc,
                         'vkeys': [u, v],
                         'end_vkeys': [u, v],
                         'int_vkeys': {}, }

    # --------------------------------------------------------------------------
    #   3.  arc intersections --> cross adjacencies
    # --------------------------------------------------------------------------
    arc_pairs_seen = set()
    for arckey_1 in arcs:
        for arckey_2 in arcs:
            if arckey_1 != arckey_2:
                arc_pair = frozenset([arckey_1, arckey_2])
                if arc_pair not in arc_pairs_seen:
                    arc_1 = arcs[arckey_1]['arc']
                    arc_2 = arcs[arckey_2]['arc']
                    intersection = _curve_curve_intx(arc_1, arc_2)
                    if intersection:

                        vkey = None
                        for key in egi.nodes():
                            if distance(egi.node_coordinates(key), intersection) < TOL.angular:
                                vkey = key
                                break

                        # if intersection is not an existing node --------------
                        if vkey is None:
                            new_vkey = max(egi.nodes()) + 1
                            new_normal = normalize_vector(subtract_vectors(intersection, origin))
                            egi.add_node(new_vkey,
                                         x=intersection[0],
                                         y=intersection[1],
                                         z=intersection[2],
                                         type='zero',
                                         normal=new_normal,
                                         magnitude=0,
                                         nbrs=[])
                            arcs[arckey_1]['vkeys'].append(new_vkey)
                            arcs[arckey_2]['vkeys'].append(new_vkey)
                            arcs[arckey_1]['int_vkeys'][new_vkey] = arckey_2
                            arcs[arckey_2]['int_vkeys'][new_vkey] = arckey_1

                        # if intersection already exists -----------------------
                        else:
                            if vkey not in arcs[arckey_1]['vkeys']:
                                arcs[arckey_1]['vkeys'].append(vkey)
                                arcs[arckey_1]['int_vkeys'][vkey] = arckey_2
                            if vkey not in arcs[arckey_2]['vkeys']:
                                arcs[arckey_2]['vkeys'].append(vkey)
                                arcs[arckey_2]['int_vkeys'][vkey] = arckey_1
                        arc_pairs_seen.add(arc_pair)

    # --------------------------------------------------------------------------
    #   5.  Reorder vertices along each arc and add edges to EGI network
    # --------------------------------------------------------------------------
    for arckey in arcs:
        vkeys = arcs[arckey]['vkeys']
        if len(vkeys) > 2:
            pt_list = [egi.node_coordinates(key) for key in vkeys]
            arcs[arckey]['vkeys'] = _reorder_pts_on_arc(pt_list,
                                                        arcs[arckey]['vkeys'],
                                                        arcs[arckey]['arc'])[1]
        for i in range(len(arcs[arckey]['vkeys']) - 1):
            vkey_1 = arcs[arckey]['vkeys'][i]
            vkey_2 = arcs[arckey]['vkeys'][i + 1]
            egi.node_attribute(vkey_1, 'nbrs', egi.node_attribute(vkey_1, 'nbrs') + [vkey_2])
            egi.node_attribute(vkey_2, 'nbrs', egi.node_attribute(vkey_2, 'nbrs') + [vkey_1])
            egi.add_edge(vkey_1, vkey_2)

    # --------------------------------------------------------------------------
    #   6.  For each vertex, sort nbrs in ccw order
    # --------------------------------------------------------------------------
    _egi_sort_v_nbrs(egi)

    # --------------------------------------------------------------------------
    #   7.  Add EGI Network faces
    # --------------------------------------------------------------------------
    egi_mesh = EGI(name='egi')
    for vkey, attr in egi.nodes(data=True):
        egi_mesh.add_vertex(key=vkey, attr_dict=attr)

    egi_mesh.attributes['origin'] = list(origin)

    _egi_find_faces(egi, egi_mesh)

    # --------------------------------------------------------------------------
    #   8.  Face normals point away from the origin
    # --------------------------------------------------------------------------
    outward = [dot_vectors(egi_mesh.face_normal(fkey), subtract_vectors(egi_mesh.face_centroid(fkey), origin)) > 0 for fkey in egi_mesh.faces()]
    if not any(outward):
        egi_mesh.flip_cycles()
    elif not all(outward):
        raise Exception('The faces of the egi are not consistently oriented.')

    return egi_mesh


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   from egi to unit polyhedron
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************

def cell_from_egi(egi: EGI) -> Cell:
    """Construct a cell from an egi.

    Parameters
    ----------
    egi
        A mesh object representing the egi.

    Returns
    -------
    Cell
        A mesh object representing the unit polyhedron.

    """
    cell = Cell(name='cell')

    for fkey in egi.faces():
        x, y, z = egi.face_center(fkey)
        cell.add_vertex(key=fkey, x=x, y=y, z=z)

    for vkey in egi.vertices():
        cell_face = egi.vertex_faces(vkey, ordered=True)
        cell.add_face(cell_face[::-1], fkey=vkey)

        cell.face_attribute(vkey, 'type', egi.vertex_attribute(vkey, 'type'))

    return cell


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   helpers
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


def _draw_arc(normal_1: Sequence[float], normal_2: Sequence[float], origin: Sequence[float]) -> tuple[list[float], list[float], list[float], list[float]]:
    n1 = normalize_vector(normal_1)
    n2 = normalize_vector(normal_2)
    return n1, n2, normalize_vector(cross_vectors(n1, n2)), list(origin)


def _reorder_pts_on_arc(pt_list: list[list[float]], pt_key_list: list[int], arc: tuple) -> tuple[list[list[float]], list[int]]:
    # all points should be on the arc...
    n1, _, _, origin = arc
    angles = [angle_vectors(n1, subtract_vectors(pt, origin)) for pt in pt_list]
    order = sorted(range(len(pt_list)), key=lambda i: angles[i])
    ordered_pt_list = [pt_list[i] for i in order]
    ordered_pt_key_list = [pt_key_list[i] for i in order]
    return ordered_pt_list, ordered_pt_key_list


def _curve_curve_intx(arc_1: tuple, arc_2: tuple) -> list[float] | None:
    n1_1, n2_1, c_1, origin = arc_1
    n1_2, n2_2, c_2, _ = arc_2
    d = cross_vectors(c_1, c_2)
    if TOL.is_zero(length_vector(d)):
        return None
    d = normalize_vector(d)
    for p in (d, scale_vector(d, -1)):
        on_arc_1 = dot_vectors(cross_vectors(n1_1, p), c_1) >= -TOL.absolute and dot_vectors(cross_vectors(p, n2_1), c_1) >= -TOL.absolute
        on_arc_2 = dot_vectors(cross_vectors(n1_2, p), c_2) >= -TOL.absolute and dot_vectors(cross_vectors(p, n2_2), c_2) >= -TOL.absolute
        if on_arc_1 and on_arc_2:
            return add_vectors(p, origin)
    return None


def _egi_sort_v_nbrs(egi: Graph) -> None:
    """ The neighbors are sorted by their angle around the vertex normal;
    only the cyclic order matters.
    """
    for vkey in egi.nodes():
        nbrs = egi.node_attribute(vkey, 'nbrs')
        point = egi.node_coordinates(vkey)
        frame = Frame.from_plane(Plane(point, egi.node_attribute(vkey, 'normal')))
        p_list = []
        for nbr_vkey in nbrs:
            vector = subtract_vectors(egi.node_coordinates(nbr_vkey), point)
            p_list.append(atan2(dot_vectors(vector, frame.yaxis), dot_vectors(vector, frame.xaxis)))
        sorted_nbrs = [key for (param, key) in sorted(zip(p_list, nbrs))]
        egi.node_attribute(vkey, 'sorted_nbrs', sorted_nbrs)


def _egi_find_edge_face(u: int, v: int, egi: Graph) -> list[int]:
    """ same as duality.algorithms.find_edge_faces... using 'sorted_nbrs' instead
    """
    cycle = [u]
    while True:
        cycle.append(v)
        nbrs = egi.node_attribute(v, 'sorted_nbrs')
        nbr = nbrs[nbrs.index(u) - 1]
        u, v = v, nbr
        if v == cycle[0]:
            cycle.append(v)
            break
    return cycle


def _egi_find_faces(egi: Graph, egi_mesh: EGI) -> EGI:
    """ Modified, and simplified version of duality.algorithms.find_network_faces...
    since there are no leaves or open faces in a egi network.
    """
    u = sorted(egi.nodes(data=True), key=lambda x: (x[1]['y'], x[1]['x']))[0][0]
    v = graph_node_find_first_neighbor(egi, u)

    halfedges = [(u, v)]
    for u, v in egi.edges():
        halfedges.append((u, v))
        halfedges.append((v, u))

    visited = set()
    for u, v in halfedges:
        if (u, v) not in visited:
            cycle = _egi_find_edge_face(u, v, egi)
            visited.update(zip(cycle[:-1], cycle[1:]))
            egi_mesh.add_face(cycle[:-1])

    return egi_mesh


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   main
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


if __name__ == '__main__':
    pass
