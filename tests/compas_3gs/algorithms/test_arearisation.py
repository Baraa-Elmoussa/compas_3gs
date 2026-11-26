import json
import os

import pytest

from compas.geometry import Box
from compas.geometry import distance_point_plane
from compas.geometry import dot_vectors
from compas.geometry import subtract_vectors
from compas.geometry import volume_polyhedron
from compas.tolerance import TOL
from compas_3gs import DATA
from compas_3gs.algorithms import cell_arearise_face
from compas_3gs.diagrams import Cell

with open(os.path.join(DATA, "cell_gss.json")) as f:
    CELL_GSS = json.load(f)

# cells[0] has vertices of degree 4; cells[1] is a square frustum: bottom 10 x 10 at z = 0 (face 4), top 2.5 x 2.5 at z = 5 (face 5),
# its sides meet in an apex at z = 6.67
with open(os.path.join(DATA, "cells.json")) as f:
    CELLS = json.load(f)["cells"]

# a square frustum with almost vertical sides: bottom 10 x 10 at z = 0, top 9.9 x 9.9 at z = 10 (face 1)
STEEP_FRUSTUM = {
    "vertices": [[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0], [0.05, 0.05, 10], [9.95, 0.05, 10], [9.95, 9.95, 10], [0.05, 9.95, 10]],
    "faces": [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]],
}

# a wedge: its bottom face (face 0, 2 x 1 at z = 0) shrinks linearly as it moves up, to nothing at z = 2, and turns over beyond
WEDGE = {
    "vertices": [[0, 0, 0], [2, 0, 0], [0, 0, 2], [0, 1, 0], [2, 1, 0], [0, 1, 2]],
    "faces": [[0, 3, 4, 1], [0, 2, 5, 3], [1, 4, 5, 2], [0, 1, 2], [3, 5, 4]],
}


@pytest.mark.parametrize(
    "data, fkey, target, zmax",
    [
        pytest.param(CELL_GSS, 0, 5.0, None, id="gss-larger"),
        pytest.param(CELL_GSS, 0, 1.0, None, id="gss-smaller"),
        pytest.param(CELL_GSS, 2, 2.6, None, id="gss-face-2"),
        pytest.param(CELLS[1], 5, 20.0, None, id="frustum-top-larger"),
        # this face shrinks by about 0.2 per unit move, so the target is beyond the default search range of 1.01
        pytest.param(STEEP_FRUSTUM, 1, 97.0, 10.0, id="steep-frustum"),
    ],
)
def test_cell_arearise_face_reaches_target_area(data, fkey, target, zmax):
    cell = Cell.from_vertices_and_faces(data["vertices"], data["faces"])

    cell_arearise_face(cell, fkey, target, tol=1e-9, zmax=zmax)

    assert abs(cell.face_area(fkey) - target) < 1e-6


def test_cell_arearise_face_tol_is_a_distance_along_the_normal():
    exact = Cell.from_vertices_and_faces(CELL_GSS["vertices"], CELL_GSS["faces"])
    cell_arearise_face(exact, 2, 2.6, tol=1e-12)
    cell = Cell.from_vertices_and_faces(CELL_GSS["vertices"], CELL_GSS["faces"])
    cell_arearise_face(cell, 2, 2.6, tol=1e-3)

    offset = dot_vectors(subtract_vectors(cell.face_centroid(2), exact.face_centroid(2)), exact.face_normal(2))
    assert abs(offset) < 1e-3


def test_cell_arearise_face_finds_the_nearer_position():
    # the top face of the frustum has an area of 2 at 0.72 and at 2.61 above its position, on either side of the apex
    cell = Cell.from_vertices_and_faces(CELLS[1]["vertices"], CELLS[1]["faces"])

    cell_arearise_face(cell, 5, 2.0, tol=1e-9)

    assert cell.face_centroid(5)[2] == pytest.approx(5 + 5 / 3 * (1 - (2 / 6.25) ** 0.5))


def test_cell_arearise_face_turns_the_face_over_for_a_negative_target():
    cell = Cell.from_vertices_and_faces(WEDGE["vertices"], WEDGE["faces"])
    normal = cell.face_normal(0)

    cell_arearise_face(cell, 0, -0.5, tol=1e-9)

    assert cell.face_area(0) == pytest.approx(0.5)
    assert dot_vectors(cell.face_normal(0), normal) < 0


def test_cell_arearise_face_raises_beyond_the_search_range():
    cell = Cell.from_vertices_and_faces(STEEP_FRUSTUM["vertices"], STEEP_FRUSTUM["faces"])

    with pytest.raises(ValueError, match="does not reach"):
        cell_arearise_face(cell, 1, 97.0)


def test_cell_arearise_face_raises_beyond_the_apex():
    # the top face of the frustum shrinks to a point at the apex and grows again beyond it, so it never turns over
    cell = Cell.from_vertices_and_faces(CELLS[1]["vertices"], CELLS[1]["faces"])

    with pytest.raises(ValueError, match="does not reach"):
        cell_arearise_face(cell, 5, -1.0, zmax=10.0)


@pytest.mark.parametrize(
    "data, fkey, target",
    [
        pytest.param(CELL_GSS, 0, 5.0, id="gss-larger"),
        pytest.param(CELL_GSS, 0, 1.0, id="gss-smaller"),
        pytest.param(CELLS[1], 5, 20.0, id="frustum-top-larger"),
        # faces with vertices of degree 4, which are split first
        pytest.param(CELLS[0], 2, 50.0, id="indeterminate-larger"),
        pytest.param(CELLS[0], 6, 15.0, id="indeterminate-smaller"),
        # two positions have the target area, on either side of the apex of the sides: the nearer one keeps the cell convex
        pytest.param(CELLS[1], 4, 50.0, id="frustum-bottom-smaller"),
        pytest.param(CELLS[1], 5, 2.0, id="frustum-top-smaller"),
    ],
)
def test_cell_arearise_face_keeps_cell_closed(data, fkey, target):
    cell = Cell.from_vertices_and_faces(data["vertices"], data["faces"])
    area = cell.face_area(fkey)
    normal = cell.face_normal(fkey)

    cell_arearise_face(cell, fkey, target)

    # closed, outward, convex like the input, with planar faces
    assert cell.is_closed()
    assert cell.number_of_vertices() - cell.number_of_edges() + cell.number_of_faces() == 2
    assert volume_polyhedron(cell.to_vertices_and_faces()) > 0
    for fkey_ in cell.faces():
        centroid = cell.face_centroid(fkey_)
        normal_ = cell.face_normal(fkey_)
        for vkey in cell.vertices():
            assert dot_vectors(subtract_vectors(cell.vertex_coordinates(vkey), centroid), normal_) < TOL.absolute
        for vkey in cell.face_vertices(fkey_):
            assert distance_point_plane(cell.vertex_coordinates(vkey), (centroid, normal_)) < TOL.absolute
    # the face moved along its normal, towards the target, and all its vertices have degree 3
    assert TOL.is_allclose(cell.face_normal(fkey), normal)
    assert (cell.face_area(fkey) - area) * (target - area) > 0
    for vkey in cell.face_vertices(fkey):
        assert cell.vertex_degree(vkey) == 3


def test_cell_arearise_face_splits_indeterminate_vertices():
    cell = Cell.from_vertices_and_faces(CELLS[0]["vertices"], CELLS[0]["faces"])
    assert max(cell.vertex_degree(vkey) for vkey in cell.face_vertices(2)) == 4
    faces = sorted(cell.faces())

    cell_arearise_face(cell, 2, 50.0)

    assert sorted(cell.faces()) == faces
    assert cell.number_of_vertices() == 16
    assert all(cell.vertex_degree(vkey) == 3 for vkey in cell.vertices())


def test_cell_arearise_face_cube_raises():
    # the faces of a box keep their area when they move along their normal
    cell = Cell.from_shape(Box(1))

    with pytest.raises(Exception, match="will not change its area"):
        cell_arearise_face(cell, 0, 2.0)


def test_cell_arearise_face_callback():
    cell = Cell.from_vertices_and_faces(CELL_GSS["vertices"], CELL_GSS["faces"])

    calls = []
    cell_arearise_face(cell, 0, 5.0, callback=lambda cell, args: calls.append(args), callback_args="args")

    assert calls
    assert all(args == "args" for args in calls)


def test_cell_arearise_face_callback_not_callable():
    cell = Cell.from_vertices_and_faces(CELL_GSS["vertices"], CELL_GSS["faces"])

    with pytest.raises(Exception, match="not callable"):
        cell_arearise_face(cell, 0, 5.0, callback=1)
