import pytest

from compas.geometry import Box
from compas.geometry import add_vectors
from compas.geometry import centroid_points
from compas.geometry import dot_vectors
from compas.geometry import normalize_vector
from compas.geometry import scale_vector
from compas.geometry import subtract_vectors
from compas.geometry import volume_polyhedron
from compas.tolerance import TOL
from compas_3gs.diagrams import Cell
from compas_3gs.operations import cell_collapse_short_edge
from compas_3gs.operations import cell_face_subdivide_barycentric
from compas_3gs.operations import cell_merge_coplanar_adjacent_faces
from compas_3gs.operations import cell_relocate_face
from compas_3gs.operations import cell_split_indet_face_vertices

# `Mesh.volume` returns the absolute value; `volume_polyhedron` keeps the sign, positive for outward faces.


@pytest.fixture
def cube():
    return Cell.from_shape(Box(1))


@pytest.fixture(params=["consecutive", "gapped"])
def octahedron(request):
    vertices, faces = Cell.from_polyhedron(8).to_vertices_and_faces()
    if request.param == "consecutive":
        return Cell.from_vertices_and_faces(vertices, faces)
    cell = Cell()
    for i, (x, y, z) in enumerate(vertices):
        cell.add_vertex(key=10 + 3 * i, x=x, y=y, z=z)
    for i, face in enumerate(faces):
        cell.add_face([10 + 3 * j for j in face], fkey=5 + 7 * i)
    return cell


# ==============================================================================
# cell vertices
# ==============================================================================


def test_cell_split_indet_face_vertices(octahedron):
    cell = octahedron
    face = list(cell.faces())[2]
    assert all(cell.vertex_degree(vertex) == 4 for vertex in cell.face_vertices(face))

    result = cell_split_indet_face_vertices(cell, face)

    assert result is cell
    # every vertex of the triangle is split in two, joined by a new zero-length edge
    assert (cell.number_of_vertices(), cell.number_of_edges(), cell.number_of_faces()) == (9, 15, 8)
    vertices = cell.face_vertices(face)
    assert len(vertices) == 6
    for vertex in vertices:
        assert cell.vertex_degree(vertex) == 3
        assert len([nbr for nbr in cell.vertex_neighbors(vertex) if nbr not in vertices]) == 1
    assert cell.is_closed()
    assert volume_polyhedron(cell.to_vertices_and_faces()) == pytest.approx(4 / 3)
    center = centroid_points(cell.vertices_attributes("xyz"))
    for fkey in cell.faces():
        assert dot_vectors(cell.face_normal(fkey), subtract_vectors(cell.face_centroid(fkey), center)) > 0


def test_cell_split_indet_face_vertices_cube(cube):
    # the vertices of a cube have degree 3: nothing to split
    faces = {fkey: cube.face_vertices(fkey) for fkey in cube.faces()}
    cell_split_indet_face_vertices(cube, 5)
    assert {fkey: cube.face_vertices(fkey) for fkey in cube.faces()} == faces
    assert cube.number_of_vertices() == 8


# ==============================================================================
# cell edges
# ==============================================================================


@pytest.fixture
def pinched():
    # a unit cube with the back edge of its top face shortened to 0.04
    cell = Cell.from_shape(Box(1))
    for vertex in cell.vertices():
        x, y, z = cell.vertex_coordinates(vertex)
        if y > 0 and z > 0:
            cell.vertex_attribute(vertex, "x", 0.02 if x > 0 else -0.02)
    return cell


def test_cell_collapse_short_edge(pinched):
    cell = pinched
    u, v = [vertex for vertex in cell.vertices() if cell.vertex_attribute(vertex, "y") > 0 and cell.vertex_attribute(vertex, "z") > 0]
    assert v in cell.vertex_neighbors(u)
    assert TOL.is_close(cell.edge_length((u, v)), 0.04)
    faces = {fkey: cell.face_vertices(fkey) for fkey in cell.faces()}

    result = cell_collapse_short_edge(cell, u, v)

    assert result is cell
    assert cell.vertex_coordinates(u) == pytest.approx([0, 0.5, 0.5])
    assert cell.vertex_coordinates(v) == pytest.approx([0, 0.5, 0.5])
    # geometric only: same vertices, same faces
    assert {fkey: cell.face_vertices(fkey) for fkey in cell.faces()} == faces
    assert cell.is_closed()
    assert volume_polyhedron(cell.to_vertices_and_faces()) > 0
    center = centroid_points(cell.vertices_attributes("xyz"))
    for fkey in cell.faces():
        assert dot_vectors(cell.face_normal(fkey), subtract_vectors(cell.face_centroid(fkey), center)) > 0


def test_cell_collapse_short_edge_fixed(pinched):
    cell = pinched
    u, v = [vertex for vertex in cell.vertices() if cell.vertex_attribute(vertex, "y") > 0 and cell.vertex_attribute(vertex, "z") > 0]
    xu = cell.vertex_attribute(u, "x")
    cell.vertex_attribute(u, "x_fix", True)

    cell_collapse_short_edge(cell, u, v)

    assert cell.vertex_coordinates(u) == pytest.approx([xu, 0.5, 0.5])
    assert cell.vertex_coordinates(v) == pytest.approx([0, 0.5, 0.5])


def test_cell_collapse_short_edge_long(cube):
    xyz = {vertex: cube.vertex_coordinates(vertex) for vertex in cube.vertices()}
    u, v = next(iter(cube.edges()))
    cell_collapse_short_edge(cube, u, v, min_length=0.5)
    assert {vertex: cube.vertex_coordinates(vertex) for vertex in cube.vertices()} == xyz


# ==============================================================================
# cell faces
# ==============================================================================


@pytest.mark.parametrize(
    "xyz, normal, volume",
    [
        ([0, 0, 1], [0, 0, 1], 1.5),
        ([0, 0, 0.25], [0, 0, 1], 0.75),
        ([0, 0, 0.75], [0.2, 0, 1], 1.25),
    ],
)
def test_cell_relocate_face_cube(cube, xyz, normal, volume):
    top = [fkey for fkey in cube.faces() if cube.face_normal(fkey)[2] > 0.9][0]

    result = cell_relocate_face(cube, top, xyz, normal)

    assert result is cube
    assert list(cube.face_centroid(top)) == pytest.approx(xyz)
    assert list(cube.face_normal(top)) == pytest.approx(normalize_vector(normal))
    for vertex in cube.face_vertices(top):
        assert TOL.is_zero(dot_vectors(subtract_vectors(cube.vertex_coordinates(vertex), xyz), normal))
    assert (cube.number_of_vertices(), cube.number_of_faces()) == (8, 6)
    assert cube.is_closed()
    assert volume_polyhedron(cube.to_vertices_and_faces()) == pytest.approx(volume)
    center = centroid_points(cube.vertices_attributes("xyz"))
    for fkey in cube.faces():
        assert dot_vectors(cube.face_normal(fkey), subtract_vectors(cube.face_centroid(fkey), center)) > 0


def test_cell_relocate_face_octahedron(octahedron):
    cell = octahedron
    face = list(cell.faces())[2]
    normal = cell.face_normal(face)
    xyz = add_vectors(cell.face_centroid(face), scale_vector(normal, 0.3))

    cell_relocate_face(cell, face, xyz, normal)

    # the face is split to a hexagon first, then every vertex slides along its one outside edge
    assert len(cell.face_vertices(face)) == 6
    assert list(cell.face_centroid(face)) == pytest.approx(xyz)
    assert list(cell.face_normal(face)) == pytest.approx(list(normal))
    for vertex in cell.face_vertices(face):
        assert TOL.is_zero(dot_vectors(subtract_vectors(cell.vertex_coordinates(vertex), xyz), normal))
    assert cell.is_closed()
    assert volume_polyhedron(cell.to_vertices_and_faces()) > 4 / 3
    center = centroid_points(cell.vertices_attributes("xyz"))
    for fkey in cell.faces():
        assert dot_vectors(cell.face_normal(fkey), subtract_vectors(cell.face_centroid(fkey), center)) > 0


def test_cell_merge_coplanar_adjacent_faces_split_cube(cube):
    for fkey in list(cube.faces()):
        vertices = cube.face_vertices(fkey)
        cube.split_face(fkey, vertices[0], vertices[2])
    assert cube.number_of_faces() == 12

    result = cell_merge_coplanar_adjacent_faces(cube)

    assert result is cube
    assert cube.number_of_faces() == 6
    assert all(len(cube.face_vertices(fkey)) == 4 for fkey in cube.faces())
    assert sorted(tuple(round(c) for c in cube.face_normal(fkey)) for fkey in cube.faces()) == sorted([(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)])
    assert sum(cube.face_area(fkey) for fkey in cube.faces()) == pytest.approx(6)
    assert cube.is_closed()
    assert volume_polyhedron(cube.to_vertices_and_faces()) == pytest.approx(1)


def test_cell_merge_coplanar_adjacent_faces_fan(cube):
    for fkey in list(cube.faces()):
        cell_face_subdivide_barycentric(cube, fkey)
    assert (cube.number_of_vertices(), cube.number_of_faces()) == (14, 24)

    cell_merge_coplanar_adjacent_faces(cube)

    # each fan of 4 triangles becomes the quad again, without its centre vertex
    assert (cube.number_of_vertices(), cube.number_of_faces()) == (8, 6)
    assert all(len(cube.face_vertices(fkey)) == 4 for fkey in cube.faces())
    assert sum(cube.face_area(fkey) for fkey in cube.faces()) == pytest.approx(6)
    assert cube.is_closed()
    assert volume_polyhedron(cube.to_vertices_and_faces()) == pytest.approx(1)


def test_cell_merge_coplanar_adjacent_faces_tol(cube):
    # split the top face and raise one corner by 0.2: the two triangles are 15.8 degrees apart
    top = [fkey for fkey in cube.faces() if cube.face_normal(fkey)[2] > 0.9][0]
    vertices = cube.face_vertices(top)
    cube.split_face(top, vertices[0], vertices[2])
    cube.vertex_attribute(vertices[1], "z", 0.7)

    cell_merge_coplanar_adjacent_faces(cube)
    assert cube.number_of_faces() == 7

    cell_merge_coplanar_adjacent_faces(cube, tol=0.1)
    assert cube.number_of_faces() == 6
    assert cube.is_closed()
    assert volume_polyhedron(cube.to_vertices_and_faces()) > 0


def test_cell_face_subdivide_barycentric(cube):
    top = [fkey for fkey in cube.faces() if cube.face_normal(fkey)[2] > 0.9][0]
    vertices = list(cube.vertices())

    result = cell_face_subdivide_barycentric(cube, top)

    assert result is cube
    assert not cube.has_face(top)
    assert (cube.number_of_vertices(), cube.number_of_faces()) == (9, 9)
    center = [vertex for vertex in cube.vertices() if vertex not in vertices]
    assert len(center) == 1
    assert cube.vertex_coordinates(center[0]) == pytest.approx([0, 0, 0.5])
    fan = cube.vertex_faces(center[0])
    assert len(fan) == 4
    for fkey in fan:
        assert len(cube.face_vertices(fkey)) == 3
        assert list(cube.face_normal(fkey)) == pytest.approx([0, 0, 1])
        assert cube.face_area(fkey) == pytest.approx(0.25)
    assert sum(cube.face_area(fkey) for fkey in cube.faces()) == pytest.approx(6)
    assert cube.is_closed()
    assert volume_polyhedron(cube.to_vertices_and_faces()) == pytest.approx(1)


def test_cell_face_subdivide_barycentric_octahedron(octahedron):
    cell = octahedron
    area = sum(cell.face_area(fkey) for fkey in cell.faces())

    for fkey in list(cell.faces()):
        cell_face_subdivide_barycentric(cell, fkey)

    assert (cell.number_of_vertices(), cell.number_of_faces()) == (14, 24)
    assert sum(cell.face_area(fkey) for fkey in cell.faces()) == pytest.approx(area)
    assert cell.is_closed()
    assert volume_polyhedron(cell.to_vertices_and_faces()) == pytest.approx(4 / 3)
    center = centroid_points(cell.vertices_attributes("xyz"))
    for fkey in cell.faces():
        assert dot_vectors(cell.face_normal(fkey), subtract_vectors(cell.face_centroid(fkey), center)) > 0
