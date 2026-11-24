
from compas.datastructures import Graph
from compas.datastructures import VolMesh

from compas.geometry import dot_vectors
from compas.geometry import normal_polygon
from compas.geometry import subtract_vectors

from compas_3gs.datastructures import Network3gs
from compas_3gs.datastructures import VolMesh3gs


__all__ = ['volmesh_dual_volmesh',
           'volmesh_dual_network']


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   volmesh  >>>  dual volmesh
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


def volmesh_dual_volmesh(volmesh: VolMesh, cls: type[VolMesh] | None = None) -> VolMesh:
    """Constructs the dual volmesh of a volmesh.

    Every cell of the volmesh becomes a vertex of the dual (vertex key = cell key), and every interior
    vertex becomes a cell of the dual (cell key = vertex key), with one face per edge of the vertex.
    The cells around an edge give the vertices of its dual face, in cyclic order around the edge,
    oriented so that the face points out of its dual cell.

    Parameters
    ----------
    volmesh
        A volmesh object.
    cls
        The class of the dual volmesh. Default is `VolMesh3gs`.

    Returns
    -------
    VolMesh
        The dual volmesh object of the input volmesh.

    Raises
    ------
    Exception
        If the volmesh has no interior vertices.

    """

    if not cls:
        cls = VolMesh3gs

    # 1. make volmesh instance -------------------------------------------------
    dual_volmesh = cls(name='volmesh_dual_volmesh')

    # 2. add vertex for each cell ----------------------------------------------
    for cell in volmesh.cells():
        x, y, z = volmesh.cell_centroid(cell)
        dual_volmesh.add_vertex(key=cell, x=x, y=y, z=z)

    # 3. find interior vertices ------------------------------------------------
    ext_vkeys = []
    boundary_hfkeys = volmesh.halffaces_on_boundaries()
    for halfface in boundary_hfkeys:
        for vertex in volmesh.halfface_vertices(halfface):
            ext_vkeys.append(vertex)
    int_vkeys = list(set(volmesh.vertices()) - set(ext_vkeys))

    if len(int_vkeys) < 1:
        raise Exception('Not enough cells to create a dual volmesh.')

    # 4. for each interior vertex, find neighbors ------------------------------
    for u in int_vkeys:
        cell_halffaces = []
        for v in volmesh.vertex_neighbors(u):

            # walk around the edge: the face of a cell with halfedge (u, v) leads to the next cell
            # (edge_cells lists the cells around the edge, but not in cyclic order)
            uv_halfface = {volmesh.halfface_cell(hfkey): hfkey for hfkey in volmesh.edge_halffaces((u, v))}
            cell = next(iter(uv_halfface))
            halfface = []
            while cell not in halfface:
                halfface.append(cell)
                cell = volmesh.halfface_opposite_cell(uv_halfface[cell])

            # the dual face points from u to v, out of the dual cell of u
            normal = normal_polygon([dual_volmesh.vertex_coordinates(ckey) for ckey in halfface])
            if dot_vectors(normal, subtract_vectors(volmesh.vertex_coordinates(v), volmesh.vertex_coordinates(u))) < 0:
                halfface.reverse()

            cell_halffaces.append(halfface)

        dual_volmesh.add_cell(cell_halffaces, ckey=u)

    return dual_volmesh


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   volmesh  >>>  dual network
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


def volmesh_dual_network(volmesh: VolMesh, cls: type[Graph] | None = None) -> Graph:
    """Computes the dual network dual of a volmesh.

    Every cell becomes a node at the cell centroid (node key = cell key), and every pair of
    neighbouring cells an edge, directed from the cell that is visited first.

    Parameters
    ----------
    volmesh
        A volmesh object.
    cls
        The class of the dual network. Default is `Network3gs`.

    Returns
    -------
    Graph
        The dual network object of the input volmesh.

    """
    if not cls:
        cls = Network3gs

    dual_network = cls()

    for cell in volmesh.cells():
        x, y, z = volmesh.cell_centroid(cell)
        dual_network.add_node(cell, x=x, y=y, z=z)

    for cell in volmesh.cells():
        # cell_neighbors lists a neighbour once per shared halfedge
        for nbr in dict.fromkeys(volmesh.cell_neighbors(cell)):
            if dual_network.has_edge((cell, nbr), directed=False):
                continue
            dual_network.add_edge(cell, nbr)

    return dual_network


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   Main
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


if __name__ == '__main__':
    pass
