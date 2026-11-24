import random

import pytest

from compas.geometry import cross_vectors
from compas.geometry import dot_vectors
from compas.geometry import subtract_vectors
from compas.geometry import volume_polyhedron
from compas.tolerance import TOL
from compas_3gs.algorithms import volmesh_dual_network
from compas_3gs.algorithms import volmesh_dual_volmesh
from compas_3gs.datastructures import Network3gs
from compas_3gs.datastructures import VolMesh3gs
from compas_3gs.diagrams import ForceVolMesh
from compas_3gs.diagrams import FormNetwork


@pytest.fixture
def meshgrid():
    """A 3x3x3 meshgrid of unit cubes: 27 cells, 64 vertices, 8 of them interior."""
    return ForceVolMesh.from_meshgrid(dx=3, nx=3)


# ==============================================================================
# volmesh_dual_volmesh
# ==============================================================================


def test_volmesh_dual_volmesh_keys_and_counts(meshgrid):
    dual = volmesh_dual_volmesh(meshgrid)

    interior = [vertex for vertex in meshgrid.vertices() if not meshgrid.is_vertex_on_boundary(vertex)]
    assert type(dual) is VolMesh3gs
    assert sorted(dual.cells()) == sorted(interior)
    assert sorted(dual.vertices()) == sorted(meshgrid.cells())
    for cell in meshgrid.cells():
        assert TOL.is_allclose(dual.vertex_coordinates(cell), meshgrid.cell_centroid(cell))

    # a 2x2x2 grid of cubes
    assert dual.number_of_cells() == 8
    assert dual.number_of_vertices() == 27
    assert dual.number_of_edges() == 54
    assert dual.number_of_faces() == 36
    assert len(list(dual.halffaces())) == 48
    assert len(dual.halffaces_on_boundaries()) == 24


@pytest.mark.parametrize("amount", [0.0, 0.2])
def test_volmesh_dual_volmesh_closed_outward_paired(meshgrid, amount):
    rng = random.Random(1)
    for vertex in sorted(meshgrid.vertices()):
        x, y, z = meshgrid.vertex_coordinates(vertex)
        meshgrid.vertex_attributes(vertex, "xyz", [x + rng.uniform(-amount, amount), y + rng.uniform(-amount, amount), z + rng.uniform(-amount, amount)])

    dual = volmesh_dual_volmesh(meshgrid)

    assert dual.number_of_cells() == 8
    for cell in dual.cells():
        mesh = dual.cell_to_mesh(cell)
        assert mesh.number_of_vertices() == 8
        assert mesh.number_of_faces() == 6
        assert mesh.is_closed()
        assert volume_polyhedron(mesh.to_vertices_and_faces()) > 0

    # 48 halffaces, of which 24 pair up inside: a halfface without an opposite would raise the boundary count
    assert len(dual.halffaces_on_boundaries()) == 24
    for halfface in dual.halffaces():
        # a proper cycle: every corner turns the same way around the face normal (a bowtie turns both ways)
        points = dual.halfface_coordinates(halfface)
        normal = dual.halfface_normal(halfface)
        for a, b, c in zip(points, points[1:] + points[:1], points[2:] + points[:2]):
            assert dot_vectors(cross_vectors(subtract_vectors(b, a), subtract_vectors(c, b)), normal) > 0

        if dual.is_halfface_on_boundary(halfface):
            continue
        opposite = dual.halfface_opposite_halfface(halfface)
        assert opposite is not None
        vertices = dual.halfface_vertices(halfface)
        reverse = dual.halfface_vertices(opposite)[::-1]
        start = vertices.index(reverse[0])
        assert vertices[start:] + vertices[:start] == reverse


def test_volmesh_dual_volmesh_faces_are_unit_squares(meshgrid):
    dual = volmesh_dual_volmesh(meshgrid)

    for halfface in dual.halffaces():
        assert len(set(dual.halfface_vertices(halfface))) == 4
        assert dual.halfface_area(halfface) == pytest.approx(1.0)
    for cell in dual.cells():
        assert volume_polyhedron(dual.cell_to_vertices_and_faces(cell)) == pytest.approx(1.0)


def test_volmesh_dual_volmesh_face_normals_follow_edges(meshgrid):
    dual = volmesh_dual_volmesh(meshgrid)

    for u in dual.cells():
        halffaces = {frozenset(dual.halfface_vertices(halfface)): halfface for halfface in dual.cell_faces(u)}
        assert len(halffaces) == 6
        for v in meshgrid.vertex_neighbors(u):
            # the dual face of edge (u, v) is made of the cells around the edge and points from u to v
            halfface = halffaces[frozenset(meshgrid.edge_cells((u, v)))]
            assert TOL.is_allclose(dual.halfface_normal(halfface), meshgrid.edge_direction((u, v)))


def test_volmesh_dual_volmesh_cls(meshgrid):
    dual = volmesh_dual_volmesh(meshgrid, cls=ForceVolMesh)

    assert type(dual) is ForceVolMesh
    assert dual.number_of_cells() == 8


def test_volmesh_dual_volmesh_without_interior_vertex_raises():
    volmesh = ForceVolMesh.from_meshgrid(dx=1, nx=1)

    with pytest.raises(Exception, match="Not enough cells"):
        volmesh_dual_volmesh(volmesh)


# ==============================================================================
# volmesh_dual_network
# ==============================================================================


def test_volmesh_dual_network(meshgrid):
    network = volmesh_dual_network(meshgrid)

    assert type(network) is Network3gs
    assert sorted(network.nodes()) == sorted(meshgrid.cells())
    for cell in meshgrid.cells():
        assert TOL.is_allclose(network.node_coordinates(cell), meshgrid.cell_centroid(cell))

    assert network.number_of_nodes() == 27
    assert network.number_of_edges() == 54

    # no edge twice, in either direction, and no loops
    edges = [frozenset(edge) for edge in network.edges()]
    assert len(set(edges)) == 54
    assert all(len(edge) == 2 for edge in edges)

    # one edge per pair of cells that share a face
    pairs = {
        frozenset((meshgrid.halfface_cell(halfface), meshgrid.halfface_opposite_cell(halfface)))
        for halfface in meshgrid.halffaces()
        if not meshgrid.is_halfface_on_boundary(halfface)
    }
    assert set(edges) == pairs

    # directed from the cell visited first; neighbouring unit cubes are one apart
    order = list(meshgrid.cells())
    for u, v in network.edges():
        assert order.index(u) < order.index(v)
        assert network.edge_length((u, v)) == pytest.approx(1.0)


def test_volmesh_dual_network_cls(meshgrid):
    network = volmesh_dual_network(meshgrid, cls=FormNetwork)

    assert type(network) is FormNetwork
    assert network.number_of_edges() == 54
