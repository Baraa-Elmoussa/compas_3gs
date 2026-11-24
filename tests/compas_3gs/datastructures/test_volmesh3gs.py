import random
from collections import Counter

import pytest

import compas
from compas.datastructures import VolMesh
from compas.geometry import dot_vectors
from compas.geometry import normalize_vector
from compas.geometry import subtract_vectors
from compas.geometry import volume_polyhedron
from compas.tolerance import TOL
from compas_3gs.datastructures import VolMesh3gs


@pytest.fixture
def grid():
    # 2 x 2 x 2 cells of size 0.5 between (0, 0, 0) and (1, 1, 1)
    return VolMesh3gs.from_meshgrid(1, 1, 1, 2, 2, 2)


@pytest.fixture
def warped():
    grid = VolMesh.from_meshgrid(1, 1, 1, 2, 2, 2)
    rng = random.Random(0)
    vertices = [[axis + rng.uniform(-0.1, 0.1) for axis in grid.vertex_coordinates(vertex)] for vertex in grid.vertices()]
    cells = [[grid.halfface_vertices(halfface) for halfface in grid.cell_faces(cell)] for cell in grid.cells()]
    return VolMesh3gs.from_vertices_and_cells(vertices, cells)


def test_constructor():
    volmesh = VolMesh3gs()
    assert isinstance(volmesh, VolMesh)
    assert volmesh.name == "VolMesh3gs"
    assert volmesh.primal is None
    assert volmesh.attributes == {"scale": 1.0}
    assert volmesh.default_face_attributes == {"_a": 0.0}
    assert volmesh.number_of_cells() == 0
    assert VolMesh3gs(name="force").name == "force"


def test_scale_is_an_attribute(grid):
    assert grid.attributes["scale"] == 1.0
    grid.attributes["scale"] = 2.5
    assert grid.attributes["scale"] == 2.5
    assert VolMesh3gs().attributes["scale"] == 1.0


@pytest.mark.parametrize("name", ["grid", "warped"])
def test_constructed_volmesh_is_valid(name, request):
    volmesh = request.getfixturevalue(name)
    assert type(volmesh) is VolMesh3gs
    assert volmesh.number_of_cells() == 8
    assert volmesh.number_of_vertices() == 27
    assert len(list(volmesh.halffaces())) == 48
    assert len(volmesh.halffaces_on_boundaries()) == 24
    for cell in volmesh.cells():
        mesh = volmesh.cell_to_mesh(cell)
        assert mesh.is_closed()
        assert volume_polyhedron(mesh.to_vertices_and_faces()) > 0
    shared = Counter(frozenset(volmesh.halfface_vertices(halfface)) for halfface in volmesh.halffaces())
    for halfface in volmesh.halffaces():
        opposite = volmesh.halfface_opposite_halfface(halfface)
        if shared[frozenset(volmesh.halfface_vertices(halfface))] == 1:
            assert opposite is None
        else:
            assert opposite is not None
            assert volmesh.halfface_opposite_halfface(opposite) == halfface


def test_vertex_normal(grid):
    for vertex in grid.vertices():
        # -1, 0 or +1 per axis for a coordinate 0, 0.5 or 1
        signs = [round(2 * axis) - 1 for axis in grid.vertex_coordinates(vertex)]
        normal = grid.vertex_normal(vertex)
        if signs == [0, 0, 0]:
            assert not grid.is_vertex_on_boundary(vertex)
            assert normal is None
        else:
            assert grid.is_vertex_on_boundary(vertex)
            assert TOL.is_close(sum(axis**2 for axis in normal), 1.0)
            assert TOL.is_allclose(normal, normalize_vector(signs))


def test_vertex_normal_points_out_of_warped_volmesh(warped):
    centroid = warped.centroid()
    for vertex in warped.vertices():
        normal = warped.vertex_normal(vertex)
        if not warped.is_vertex_on_boundary(vertex):
            assert normal is None
            continue
        assert TOL.is_close(sum(axis**2 for axis in normal), 1.0)
        assert dot_vectors(normal, subtract_vectors(warped.vertex_coordinates(vertex), centroid)) > 0


@pytest.mark.parametrize("axis", [0, 1, 2])
def test_vertex_update_xyz_respects_fix(grid, axis):
    grid.vertex_attribute(26, "xyz"[axis] + "_fix", True)

    grid.vertex_update_xyz(26, [5.0, 6.0, 7.0])
    expected = [5.0, 6.0, 7.0]
    expected[axis] = 1.0
    assert TOL.is_allclose(grid.vertex_coordinates(26), expected)

    # a vertex without fix attributes is free
    grid.vertex_update_xyz(0, [5.0, 6.0, 7.0])
    assert TOL.is_allclose(grid.vertex_coordinates(0), [5.0, 6.0, 7.0])

    grid.vertex_update_xyz(26, [8.0, 9.0, 10.0], constrained=False)
    assert TOL.is_allclose(grid.vertex_coordinates(26), [8.0, 9.0, 10.0])


def test_cell_pair_halffaces_on_every_neighbouring_pair(warped):
    pairs = 0
    for cell in warped.cells():
        for nbr in dict.fromkeys(warped.cell_neighbors(cell)):
            halfface, opposite = warped.cell_pair_halffaces(cell, nbr)
            pairs += 1
            assert halfface in warped.cell_faces(cell)
            assert opposite in warped.cell_faces(nbr)
            assert warped.halfface_cell(halfface) == cell
            assert warped.halfface_cell(opposite) == nbr
            assert warped.halfface_opposite_halfface(halfface) == opposite
            assert warped.cell_pair_halffaces(nbr, cell) == (opposite, halfface)
            # the same face, run the other way round
            cycle = warped.halfface_vertices(halfface)
            other = warped.halfface_vertices(opposite)
            start = other.index(cycle[0])
            assert other[start:] + other[:start] == [cycle[0]] + cycle[:0:-1]
            assert TOL.is_allclose(warped.face_normal(opposite), warped.face_normal(halfface).scaled(-1))
    # 12 interior faces, each found from both sides
    assert pairs == 24


@pytest.mark.parametrize("pair", [(0, 0), (0, 3), (0, 7), (1, 2), (3, 4)])
def test_cell_pair_halffaces_raises_without_shared_face(grid, pair):
    # a cell and itself, and cells that share only an edge or a vertex
    with pytest.raises(KeyError):
        grid.cell_pair_halffaces(*pair)


@pytest.mark.roundtrip
@pytest.mark.xfail(
    strict=True,
    reason="compas 2.15.1 VolMesh.__data__ writes each face of a cell once per halfedge: a 2x2x2 meshgrid comes back with 192 halffaces instead of 48",
)
@pytest.mark.parametrize("how", ["json", "copy"])
def test_roundtrip(grid, how):
    grid.name = "grid"
    grid.attributes["scale"] = 2.0
    grid.vertex_attribute(0, "x_fix", True)
    grid.face_attribute(0, "_a", 5.0)

    if how == "json":
        other = compas.json_loads(compas.json_dumps(grid))
    else:
        other = grid.copy(copy_guid=True)

    assert type(other) is VolMesh3gs
    assert other.guid == grid.guid
    assert other.name == "grid"
    assert other.attributes == {"scale": 2.0}
    assert other.default_vertex_attributes == grid.default_vertex_attributes
    assert other.default_face_attributes == grid.default_face_attributes
    assert other.default_cell_attributes == grid.default_cell_attributes
    assert list(other.vertices()) == list(grid.vertices())
    assert list(other.cells()) == list(grid.cells())
    assert list(other.halffaces()) == list(grid.halffaces())
    for vertex in grid.vertices():
        assert TOL.is_allclose(other.vertex_coordinates(vertex), grid.vertex_coordinates(vertex))
        assert other.vertex_attributes(vertex) == grid.vertex_attributes(vertex)
    for cell in grid.cells():
        assert other.cell_faces(cell) == grid.cell_faces(cell)
    for halfface in grid.halffaces():
        assert other.halfface_vertices(halfface) == grid.halfface_vertices(halfface)
        assert other.face_attributes(halfface) == grid.face_attributes(halfface)
