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

# the default convergence tolerance of cell_arearise_face, documented as an area
TOL_AREA = 1e-3


@pytest.mark.parametrize(
    "data, fkey, target",
    [
        pytest.param(CELL_GSS, 0, 5.0, id="gss-larger"),
        pytest.param(CELL_GSS, 0, 1.0, id="gss-smaller"),
        pytest.param(CELLS[1], 5, 20.0, id="frustum-top-larger"),
        pytest.param(
            CELL_GSS,
            2,
            2.6,
            id="gss-face-2",
            marks=pytest.mark.xfail(
                strict=True,
                reason=(
                    "arearisation.py:87: `tol` (documented at arearisation.py:31 as an area) is passed to golden_section_search, "
                    "where it is a width of the move distance; face 2 of cell_gss.json goes 2.527 -> 2.5956, an area error of 0.0044 > 1e-3"
                ),
            ),
        ),
        pytest.param(
            STEEP_FRUSTUM,
            1,
            97.0,
            id="steep-frustum",
            marks=pytest.mark.xfail(
                strict=True,
                reason=(
                    "arearisation.py:84-85: the search interval [0, |area - target|] uses an area as a distance; this face shrinks by about 0.2 "
                    "per unit move, so the target is out of reach: 98.01 -> 97.81 instead of 97"
                ),
            ),
        ),
    ],
)
def test_cell_arearise_face_reaches_target_area(data, fkey, target):
    cell = Cell.from_vertices_and_faces(data["vertices"], data["faces"])

    cell_arearise_face(cell, fkey, target)

    assert abs(cell.face_area(fkey) - target) < TOL_AREA


@pytest.mark.parametrize(
    "data, fkey, target",
    [
        pytest.param(CELL_GSS, 0, 5.0, id="gss-larger"),
        pytest.param(CELL_GSS, 0, 1.0, id="gss-smaller"),
        pytest.param(CELLS[1], 5, 20.0, id="frustum-top-larger"),
        # faces with vertices of degree 4, which are split first
        pytest.param(CELLS[0], 2, 50.0, id="indeterminate-larger"),
        pytest.param(CELLS[0], 6, 15.0, id="indeterminate-smaller"),
        pytest.param(
            CELLS[1],
            4,
            50.0,
            id="frustum-bottom-smaller",
            marks=pytest.mark.xfail(
                strict=True,
                reason=(
                    "arearisation.py:84-87: the search interval [0, |area - target|] = [0, 50] reaches past the apex of the sides (a move of 6.67), "
                    "where the area grows again, and the search returns that solution: the bottom face ends at z = 11.38 instead of 1.95, "
                    "above the top face, and the cell is inverted (signed volume -82.06)"
                ),
            ),
        ),
        pytest.param(
            CELLS[1],
            5,
            2.0,
            id="frustum-top-smaller",
            marks=pytest.mark.xfail(
                strict=True,
                reason=(
                    "arearisation.py:84-87: the search interval [0, |area - target|] = [0, 4.25] reaches past the apex of the sides (a move of 1.67), "
                    "where the area grows again, and the search returns that solution: the top face ends at z = 7.61 instead of 5.72 "
                    "and the cell intersects itself (a vertex 1.13 in front of a face plane)"
                ),
            ),
        ),
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
