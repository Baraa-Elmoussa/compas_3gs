import pytest

from compas.datastructures import Mesh
from compas.geometry import Box
from compas.geometry import area_polygon
from compas.geometry import centroid_points
from compas.geometry import distance_point_point
from compas.tolerance import TOL
from compas_3gs.datastructures import VolMesh3gs
from compas_3gs.diagrams import Cell
from compas_3gs.utilities import cell_face_areaness
from compas_3gs.utilities import cell_face_flatness
from compas_3gs.utilities import datastructure_centroid
from compas_3gs.utilities import polygon_area_footprint
from compas_3gs.utilities import polygon_area_oriented
from compas_3gs.utilities import polygon_normal_oriented
from compas_3gs.utilities import resultant_vector
from compas_3gs.utilities import scale_polygon
from compas_3gs.utilities import volmesh_face_areaness
from compas_3gs.utilities import volmesh_face_flatness

SQUARE = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]]
L_SHAPE = [[0, 0, 0], [2, 0, 0], [2, 1, 0], [1, 1, 0], [1, 2, 0], [0, 2, 0]]
BOWTIE = [[0, 0, 0], [1, 1, 0], [1, 0, 0], [0, 1, 0]]
SADDLE = [[0, 0, 0], [1, 0, 0.75], [1, 1, 0], [0, 1, 0.75]]


# ==============================================================================
# vectors
# ==============================================================================


def test_resultant_vector():
    vectors = {0: [1, 0, 0], 1: [0, 2, 0], 2: [0, 0, -3]}
    locations = {0: [0, 0, 0], 1: [1, 0, 0], 2: [0, 1, 0]}
    resultant, _ = resultant_vector(vectors, locations)
    assert resultant == pytest.approx([1, 2, -3])


def test_resultant_vector_location():
    # two parallel forces of 1 and 3 at x = 0 and x = 4: the resultant acts at x = 3
    vectors = {0: [0, 0, 1], 1: [0, 0, 3]}
    locations = {0: [0, 0, 0], 1: [4, 0, 0]}
    resultant, location = resultant_vector(vectors, locations)
    assert resultant == pytest.approx([0, 0, 4])
    assert location == pytest.approx([3, 0, 0])


# ==============================================================================
# polygons
# ==============================================================================


@pytest.mark.parametrize(
    "polygon, normal, area",
    [
        (SQUARE, [0, 0, 1], 1.0),
        (SQUARE[::-1], [0, 0, -1], 1.0),
        ([[0, 0, 0], [0, 1, 0], [0, 1, 1], [0, 0, 1]], [1, 0, 0], 1.0),
        (L_SHAPE, [0, 0, 1], 3.0),
        (L_SHAPE[::-1], [0, 0, -1], 3.0),
    ],
)
def test_polygon_normal_oriented(polygon, normal, area):
    assert polygon_normal_oriented(polygon) == pytest.approx(normal)
    assert polygon_normal_oriented(polygon, unitized=False) == pytest.approx([area * n for n in normal])
    assert polygon_area_oriented(polygon) == pytest.approx(area)


def test_polygon_area_oriented_warped():
    # the vector area of a warped quad is the area of its projection on the best-fit plane
    assert polygon_normal_oriented(SADDLE) == pytest.approx([0, 0, 1])
    assert polygon_area_oriented(SADDLE) == pytest.approx(1.0)


def test_polygon_area_oriented_bowtie():
    # the two lobes of a bowtie have opposite orientation and cancel
    assert TOL.is_zero(polygon_area_oriented(BOWTIE))


@pytest.mark.parametrize(
    "polygon, area",
    [
        (SQUARE, 1.0),
        (L_SHAPE, 3.0),
        # four triangles to the centroid, each of area sqrt(1 + 0.75 ** 2) / 4
        (SADDLE, 1.25),
        # the two lobes count, each of area 1/4
        (BOWTIE, 0.5),
    ],
)
def test_polygon_area_footprint(polygon, area):
    assert polygon_area_footprint(polygon) == pytest.approx(area)


@pytest.mark.parametrize("scale", [2.0, 0.5, 0.0, -1.0])
def test_scale_polygon(scale):
    triangle = [[0, 0, 0], [3, 0, 0], [0, 3, 0]]
    center = centroid_points(triangle)
    scaled = scale_polygon(triangle, scale)

    assert len(scaled) == 3
    assert centroid_points(scaled) == pytest.approx(center)
    for point, new in zip(triangle, scaled):
        assert distance_point_point(center, new) == pytest.approx(abs(scale) * distance_point_point(center, point))
    assert area_polygon(scaled) == pytest.approx(scale**2 * 4.5)
    if scale > 0:
        assert polygon_normal_oriented(scaled) == pytest.approx([0, 0, 1])


# ==============================================================================
# mesh
# ==============================================================================


def test_cell_face_flatness():
    cell = Cell.from_shape(Box(1))
    assert cell_face_flatness(cell) == {face: pytest.approx(0.0) for face in cell.faces()}

    # raise two opposite corners of the top face by 0.2: a saddle at 0.1 from its best-fit plane z = 0.6;
    # the side faces stay in their planes
    top = [face for face in cell.faces() if cell.face_normal(face)[2] > 0.9][0]
    vertices = cell.face_vertices(top)
    for vertex in vertices[1::2]:
        cell.vertex_attribute(vertex, "z", 0.7)

    flatness = cell_face_flatness(cell)
    assert set(flatness) == set(cell.faces())
    for face, value in flatness.items():
        assert value == pytest.approx(0.1 if face == top else 0.0, abs=TOL.absolute)


def test_cell_face_areaness():
    cell = Cell.from_shape(Box(1))
    faces = list(cell.faces())
    areaness = cell_face_areaness(cell, {faces[0]: 2.0, faces[1]: 0.25})
    assert areaness == {faces[0]: pytest.approx(1.0), faces[1]: pytest.approx(0.75)}


# ==============================================================================
# volmesh
# ==============================================================================


def test_volmesh_face_flatness():
    volmesh = VolMesh3gs.from_meshgrid(2, 1, 1, 2, 1, 1)
    flatness = volmesh_face_flatness(volmesh)
    # one entry per face: 2 x 6 halffaces, the interior face once
    assert len(flatness) == 11
    assert set(flatness) == set(volmesh.faces())
    assert all(TOL.is_zero(value) for value in flatness.values())

    # push two opposite corners of the interior face x = 1 by 0.2 along x: a saddle at 0.1 from x = 1.1;
    # every other face keeps its plane
    for vertex in volmesh.vertices():
        if volmesh.vertex_coordinates(vertex) in ([1, 0, 0], [1, 1, 1]):
            volmesh.vertex_attribute(vertex, "x", 1.2)

    flatness = volmesh_face_flatness(volmesh)
    warped = [face for face, value in flatness.items() if not TOL.is_zero(value)]
    assert len(warped) == 1
    assert not volmesh.is_halfface_on_boundary(warped[0])
    assert flatness[warped[0]] == pytest.approx(0.1)


def test_volmesh_face_areaness():
    volmesh = VolMesh3gs.from_meshgrid(1, 1, 1, 2, 2, 2)
    interior = [hf for hf in volmesh.halffaces() if not volmesh.is_halfface_on_boundary(hf)][0]
    opposite = volmesh.halfface_opposite_halfface(interior)
    boundary = volmesh.halffaces_on_boundaries()[0]

    areaness = volmesh_face_areaness(volmesh, {interior: 0.5, opposite: 0.5, boundary: 0.2})
    assert areaness == {interior: pytest.approx(0.25), opposite: pytest.approx(0.25), boundary: pytest.approx(0.05)}


# ==============================================================================
# datastructures
# ==============================================================================


def test_datastructure_centroid_mesh():
    # a unit square and a 2 x 1 rectangle side by side: the vertex average is not the area centroid
    mesh = Mesh.from_vertices_and_faces([[0, 0, 0], [1, 0, 0], [3, 0, 0], [3, 1, 0], [1, 1, 0], [0, 1, 0]], [[0, 1, 4, 5], [1, 2, 3, 4]])
    assert datastructure_centroid(mesh) == pytest.approx([4 / 3, 0.5, 0])


def test_datastructure_centroid_volmesh():
    volmesh = VolMesh3gs.from_meshgrid(1, 1, 1, 2, 2, 2)
    assert datastructure_centroid(volmesh) == pytest.approx([0.5, 0.5, 0.5])

    volmesh.vertex_attribute(0, "x", -2.7)
    assert datastructure_centroid(volmesh) == pytest.approx([0.4, 0.5, 0.5])
