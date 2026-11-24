import random

import pytest

import compas
from compas.datastructures import VolMesh
from compas.geometry import Point
from compas.geometry import add_vectors
from compas.geometry import distance_point_point
from compas.geometry import dot_vectors
from compas.geometry import subtract_vectors
from compas.geometry import volume_polyhedron
from compas.tolerance import TOL
from compas_3gs.datastructures import Mesh3gs
from compas_3gs.diagrams import EGI


@pytest.fixture
def warped():
    # a 2 x 2 x 2 meshgrid with every vertex moved by up to 0.1 (cells of size 0.5)
    volmesh = VolMesh.from_meshgrid(1, 1, 1, 2, 2, 2)
    rng = random.Random(0)
    for vertex in volmesh.vertices():
        xyz = [axis + rng.uniform(-0.1, 0.1) for axis in volmesh.vertex_coordinates(vertex)]
        volmesh.vertex_attributes(vertex, "xyz", xyz)
    return volmesh


def test_constructor():
    egi = EGI()
    assert isinstance(egi, Mesh3gs)
    assert egi.name == "EGI"
    assert egi.attributes == {"origin": (0, 0, 0)}
    assert egi.default_vertex_attributes == {
        "x": 0.0,
        "y": 0.0,
        "z": 0.0,
        "x_fix": False,
        "y_fix": False,
        "z_fix": False,
        "type": None,
        "normal": None,
        "target_area": None,
    }
    assert EGI(name="e").name == "e"


@pytest.mark.parametrize("cell", range(8))
def test_from_volmesh_cell_is_an_outward_octahedron(warped, cell):
    egi = EGI.from_volmesh_cell(cell, warped)
    origin = egi.attributes["origin"]
    halffaces = warped.cell_faces(cell)

    assert type(egi) is EGI
    assert egi.name == str(cell)
    assert TOL.is_allclose(origin, warped.cell_centroid(cell))

    # one vertex per halfface, one face per cell vertex, keyed alike
    assert egi.number_of_vertices() == 6
    assert egi.number_of_faces() == 8
    assert egi.number_of_edges() == 12
    assert sorted(egi.vertices()) == sorted(halffaces)
    assert sorted(egi.faces()) == sorted(warped.cell_vertices(cell))

    # every vertex is the unit normal of its halfface, placed at the origin
    for halfface in egi.vertices():
        xyz = egi.vertex_coordinates(halfface)
        assert TOL.is_close(distance_point_point(origin, xyz), 1.0)
        assert TOL.is_allclose(xyz, add_vectors(origin, warped.face_normal(halfface)))

    # every face is the triangle of the halffaces around its cell vertex
    for vertex in egi.faces():
        assert len(egi.face_vertices(vertex)) == 3
        assert set(egi.face_vertices(vertex)) == {halfface for halfface in halffaces if vertex in warped.halfface_vertices(halfface)}

    assert egi.is_closed()
    assert egi.is_manifold()
    assert volume_polyhedron(egi.to_vertices_and_faces()) > 0
    for face in egi.faces():
        assert dot_vectors(egi.face_normal(face), subtract_vectors(egi.face_centroid(face), origin)) > 0


@pytest.mark.roundtrip
@pytest.mark.parametrize("how", ["json", "copy"])
def test_roundtrip(warped, how):
    egi = EGI.from_volmesh_cell(7, warped)  # vertex keys 42 to 47, face keys are cell vertex keys
    egi.vertex_attribute(42, "target_area", 0.25)
    egi.vertex_attribute(43, "type", "zero")

    if how == "json":
        other = compas.json_loads(compas.json_dumps(egi))
    else:
        other = egi.copy(copy_guid=True)

    assert type(other) is EGI
    assert other.guid == egi.guid
    assert other.name == "7"
    assert isinstance(other.attributes["origin"], Point)
    assert TOL.is_allclose(other.attributes["origin"], egi.attributes["origin"])
    assert other.attributes.keys() == egi.attributes.keys()
    assert other.default_vertex_attributes == egi.default_vertex_attributes
    assert other.default_edge_attributes == egi.default_edge_attributes
    assert other.default_face_attributes == egi.default_face_attributes
    assert list(other.vertices()) == list(egi.vertices())
    assert sorted(other.vertices()) == [42, 43, 44, 45, 46, 47]
    assert list(other.faces()) == list(egi.faces())
    for vertex in egi.vertices():
        assert TOL.is_allclose(other.vertex_coordinates(vertex), egi.vertex_coordinates(vertex))
        assert other.vertex_attributes(vertex) == egi.vertex_attributes(vertex)
    for face in egi.faces():
        assert other.face_vertices(face) == egi.face_vertices(face)
    assert other.is_closed()
    assert volume_polyhedron(other.to_vertices_and_faces()) > 0
