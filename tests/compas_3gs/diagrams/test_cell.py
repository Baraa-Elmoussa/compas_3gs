import pytest

import compas
from compas.geometry import Box
from compas.geometry import Frame
from compas.geometry import centroid_points
from compas.geometry import volume_polyhedron
from compas.tolerance import TOL
from compas_3gs.datastructures import Mesh3gs
from compas_3gs.diagrams import Cell


def test_constructor():
    cell = Cell()
    assert isinstance(cell, Mesh3gs)
    assert cell.name == "Cell"
    assert cell.attributes == {}
    assert cell.default_vertex_attributes == {"x": 0.0, "y": 0.0, "z": 0.0, "x_fix": False, "y_fix": False, "z_fix": False}
    assert cell.default_edge_attributes == {"target_vector": None, "target_length": None}
    assert cell.default_face_attributes == {"target_area": None, "target_normal": None}
    assert Cell(name="c").name == "c"


def test_from_vertices_and_faces():
    vertices, faces = Box(1).to_vertices_and_faces()
    cell = Cell.from_vertices_and_faces(vertices, faces)
    assert type(cell) is Cell
    assert cell.is_closed()
    assert TOL.is_close(volume_polyhedron(cell.to_vertices_and_faces()), 1.0)
    assert cell.vertex_attribute(0, "x_fix") is False
    assert cell.face_attribute(0, "target_area") is None
    assert cell.edge_attribute((0, 1), "target_length") is None


def test_cell_center_of_unit_cube_with_gapped_keys():
    vertices, faces = Box(1).to_vertices_and_faces()
    cell = Cell()
    for index, (x, y, z) in enumerate(vertices):
        cell.add_vertex(key=10 + 3 * index, x=x, y=y, z=z)
    for index, face in enumerate(faces):
        cell.add_face([10 + 3 * vertex for vertex in face], fkey=100 + 7 * index)
    assert list(cell.vertices()) == [10, 13, 16, 19, 22, 25, 28, 31]
    assert list(cell.faces()) == [100, 107, 114, 121, 128, 135]
    assert cell.is_closed()
    assert TOL.is_allclose(cell.cell_center(), [0.0, 0.0, 0.0])


def test_cell_center_is_the_volume_centroid():
    cell = Cell.from_shape(Box(1, 2, 3, frame=Frame([1, 2, 3])))
    top = next(face for face in cell.faces() if TOL.is_close(cell.face_normal(face)[2], 1.0))
    cell.insert_vertex(top, xyz=cell.face_centroid(top))
    assert cell.is_closed()
    # the extra vertex on the top face moves the vertex average, not the volume centroid
    assert centroid_points(cell.vertices_attributes("xyz"))[2] > 3.1
    assert TOL.is_allclose(cell.cell_center(), [1.0, 2.0, 3.0])


@pytest.mark.roundtrip
@pytest.mark.parametrize("how", ["json", "copy"])
def test_roundtrip(how):
    cell = Cell.from_shape(Box(1))
    cell.insert_vertex(0)  # face 0 becomes four triangles: face keys 1 to 9
    cell.name = "split cube"
    cell.attributes["note"] = "test"
    cell.vertex_attribute(0, "z_fix", True)
    cell.edge_attribute((1, 2), "target_length", 2.0)
    cell.edge_attribute((2, 3), "target_vector", [0.0, 1.0, 0.0])
    cell.face_attribute(1, "target_area", 0.5)
    cell.face_attribute(2, "target_normal", [0.0, -1.0, 0.0])

    if how == "json":
        other = compas.json_loads(compas.json_dumps(cell))
    else:
        other = cell.copy(copy_guid=True)

    assert type(other) is Cell
    assert other.guid == cell.guid
    assert other.name == "split cube"
    assert other.attributes == {"note": "test"}
    assert other.default_vertex_attributes == cell.default_vertex_attributes
    assert other.default_edge_attributes == cell.default_edge_attributes
    assert other.default_face_attributes == cell.default_face_attributes
    assert list(other.vertices()) == list(cell.vertices())
    assert list(other.faces()) == list(range(1, 10))
    assert list(other.faces()) == list(cell.faces())
    for vertex in cell.vertices():
        assert TOL.is_allclose(other.vertex_coordinates(vertex), cell.vertex_coordinates(vertex))
        assert other.vertex_attributes(vertex) == cell.vertex_attributes(vertex)
    for face in cell.faces():
        assert other.face_vertices(face) == cell.face_vertices(face)
        assert other.face_attributes(face) == cell.face_attributes(face)
    for edge in cell.edges():
        assert other.edge_attributes(edge) == cell.edge_attributes(edge)
    assert other.is_closed()
    assert TOL.is_allclose(other.cell_center(), [0.0, 0.0, 0.0])
