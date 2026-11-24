import pytest

import compas
from compas.datastructures import Mesh
from compas.geometry import Box
from compas.geometry import Vector
from compas.geometry import area_triangle
from compas.geometry import is_coplanar
from compas.geometry import volume_polyhedron
from compas.tolerance import TOL
from compas_3gs.datastructures import Mesh3gs


def test_constructor():
    mesh = Mesh3gs()
    assert isinstance(mesh, Mesh)
    assert mesh.name == "Mesh3gs"
    assert mesh.attributes == {}
    assert mesh.default_vertex_attributes == {"x": 0.0, "y": 0.0, "z": 0.0}
    assert mesh.number_of_vertices() == 0
    assert mesh.number_of_faces() == 0
    assert Mesh3gs(name="egi").name == "egi"


def test_cube_has_outward_normals_and_unit_areas():
    mesh = Mesh3gs.from_shape(Box(1))
    assert type(mesh) is Mesh3gs
    assert mesh.is_closed()
    assert mesh.number_of_faces() == 6
    assert TOL.is_close(volume_polyhedron(mesh.to_vertices_and_faces()), 1.0)
    for face in mesh.faces():
        normal = mesh.face_normal(face)
        assert isinstance(normal, Vector)
        assert TOL.is_close(normal.length, 1.0)
        # the face centroid of a unit cube at the origin is half the outward normal
        assert TOL.is_allclose(normal, [2 * axis for axis in mesh.face_centroid(face)])
        assert TOL.is_close(mesh.face_area(face), 1.0)
        assert TOL.is_allclose(mesh.face_normal(face, unitized=False), normal)


@pytest.mark.parametrize("height", [0.3, 1.0, 5.0])
def test_face_area_is_vector_area_on_warped_quad(height):
    mesh = Mesh3gs.from_vertices_and_faces([[0, 0, 0], [1, 0, height], [1, 1, 0], [0, 1, height]], [[0, 1, 2, 3]])
    points = mesh.face_coordinates(0)
    assert not is_coplanar(points)
    # the diagonals are (1, 1, 0) and (-1, 1, 0) for every height: vector area 1, normal +z,
    # while the two triangles of the quad have a larger total area
    assert TOL.is_close(mesh.face_area(0), 1.0)
    assert area_triangle(points[:3]) + area_triangle([points[0], points[2], points[3]]) > 1.0
    assert TOL.is_allclose(mesh.face_normal(0), [0, 0, 1])
    assert TOL.is_allclose(mesh.face_normal(0, unitized=False), [0, 0, 1])


@pytest.mark.parametrize("unitized", [True, False])
def test_face_normal_fallback_on_bowtie(unitized):
    mesh = Mesh3gs.from_vertices_and_faces([[0, 0, 0], [1, 1, 0], [1, 0, 0], [0, 1, 0]], [[0, 1, 2, 3]])
    assert TOL.is_zero(mesh.face_area(0))
    normal = mesh.face_normal(0, unitized=unitized)
    assert isinstance(normal, Vector)
    # (1, 1, 0) x (0, -1, 0), unitized whatever `unitized` is
    assert TOL.is_allclose(normal, [0, 0, -1])


@pytest.mark.parametrize("axis", [0, 1, 2])
def test_vertex_update_xyz_respects_fix(axis):
    mesh = Mesh3gs.from_shape(Box(1))
    mesh.vertex_attribute(0, "xyz"[axis] + "_fix", True)
    start = mesh.vertex_coordinates(0)

    mesh.vertex_update_xyz(0, [5.0, 6.0, 7.0])
    expected = [5.0, 6.0, 7.0]
    expected[axis] = start[axis]
    assert TOL.is_allclose(mesh.vertex_coordinates(0), expected)

    # a vertex without fix attributes is free
    mesh.vertex_update_xyz(1, [5.0, 6.0, 7.0])
    assert TOL.is_allclose(mesh.vertex_coordinates(1), [5.0, 6.0, 7.0])

    mesh.vertex_update_xyz(0, [8.0, 9.0, 10.0], constrained=False)
    assert TOL.is_allclose(mesh.vertex_coordinates(0), [8.0, 9.0, 10.0])


@pytest.mark.roundtrip
@pytest.mark.parametrize("how", ["json", "copy"])
def test_roundtrip(how):
    mesh = Mesh3gs.from_shape(Box(1))
    mesh.insert_vertex(0)  # face 0 becomes four triangles: face keys 1 to 9
    mesh.name = "split cube"
    mesh.attributes["note"] = "test"
    mesh.vertex_attribute(0, "x_fix", True)
    mesh.edge_attribute((1, 2), "target_length", 2.0)
    mesh.face_attribute(1, "target_area", 0.5)

    if how == "json":
        other = compas.json_loads(compas.json_dumps(mesh))
    else:
        other = mesh.copy(copy_guid=True)

    assert type(other) is Mesh3gs
    assert other.guid == mesh.guid
    assert other.name == "split cube"
    assert other.attributes == {"note": "test"}
    assert other.default_vertex_attributes == mesh.default_vertex_attributes
    assert other.default_edge_attributes == mesh.default_edge_attributes
    assert other.default_face_attributes == mesh.default_face_attributes
    assert list(other.vertices()) == list(mesh.vertices())
    assert list(other.faces()) == list(range(1, 10))
    assert list(other.faces()) == list(mesh.faces())
    for vertex in mesh.vertices():
        assert TOL.is_allclose(other.vertex_coordinates(vertex), mesh.vertex_coordinates(vertex))
        assert other.vertex_attributes(vertex) == mesh.vertex_attributes(vertex)
    for face in mesh.faces():
        assert other.face_vertices(face) == mesh.face_vertices(face)
        assert other.face_attributes(face) == mesh.face_attributes(face)
    for edge in mesh.edges():
        assert other.edge_attributes(edge) == mesh.edge_attributes(edge)
    assert other.is_closed()
    assert TOL.is_close(volume_polyhedron(other.to_vertices_and_faces()), 1.0)
