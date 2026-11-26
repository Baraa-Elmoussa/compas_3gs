import json
import os
import random
from itertools import combinations
from math import cos
from math import pi
from math import sin

import pytest

from compas.geometry import add_vectors
from compas.geometry import angle_vectors
from compas.geometry import convex_hull
from compas.geometry import cross_vectors
from compas.geometry import distance_point_plane
from compas.geometry import distance_point_point
from compas.geometry import dot_vectors
from compas.geometry import normalize_vector
from compas.geometry import subtract_vectors
from compas.geometry import volume_polyhedron
from compas.tolerance import TOL
from compas_3gs import DATA
from compas_3gs.algorithms import cell_from_egi
from compas_3gs.algorithms import egi_from_vectors
from compas_3gs.diagrams import EGI
from compas_3gs.diagrams import Cell

AXES = {0: [1, 0, 0], 1: [-1, 0, 0], 2: [0, 1, 0], 3: [0, -1, 0], 4: [0, 0, 1], 5: [0, 0, -1]}
TETRAHEDRON = {0: [1, 1, 1], 1: [1, -1, -1], 2: [-1, 1, -1], 3: [-1, -1, 1]}

with open(os.path.join(DATA, "equilibrated_force_vectors.json")) as f:
    FORCE_GROUPS = json.load(f)


@pytest.fixture
def random_vectors(request):
    """`request.param` unit vectors from a fixed seed, at least 0.05 rad apart, with the origin strictly inside their convex hull."""
    n = request.param
    rng = random.Random(n)
    while True:
        vectors = [normalize_vector([rng.gauss(0, 1) for _ in range(3)]) for _ in range(n)]
        if min(angle_vectors(u, v) for u, v in combinations(vectors, 2)) < 0.05:
            continue
        # the faces of `convex_hull` point outwards: the origin is inside if it lies behind every face, here by at least 0.01
        hull = convex_hull(vectors)
        normals = [normalize_vector(cross_vectors(subtract_vectors(vectors[b], vectors[a]), subtract_vectors(vectors[c], vectors[a]))) for a, b, c in hull]
        if all(dot_vectors(normal, vectors[face[0]]) > 0.01 for normal, face in zip(normals, hull)):
            return dict(enumerate(vectors))


# ==============================================================================
# egi_from_vectors
# ==============================================================================


def test_egi_from_vectors_octahedron():
    # the lengths of the vectors do not matter, only their directions
    vectors = {key: [2.5 * x, 2.5 * y, 2.5 * z] if key % 2 else [x, y, z] for key, (x, y, z) in AXES.items()}
    origin = [1.0, -2.0, 0.5]

    egi = egi_from_vectors(vectors, origin)

    assert type(egi) is EGI
    assert egi.attributes["origin"] == origin
    assert egi.number_of_vertices() == 6
    assert egi.number_of_edges() == 12
    assert egi.number_of_faces() == 8
    assert egi.is_closed()
    assert sorted(egi.vertices()) == sorted(AXES)
    for vkey in egi.vertices():
        assert egi.vertex_attribute(vkey, "type") == "face"
        assert TOL.is_allclose(egi.vertex_attribute(vkey, "normal"), AXES[vkey])
        assert TOL.is_allclose(egi.vertex_coordinates(vkey), add_vectors(origin, AXES[vkey]))
        assert egi.vertex_degree(vkey) == 4
    for fkey in egi.faces():
        assert len(egi.face_vertices(fkey)) == 3
        assert dot_vectors(egi.face_normal(fkey), subtract_vectors(egi.face_centroid(fkey), origin)) > 0
    # the regular octahedron inscribed in the unit sphere
    assert volume_polyhedron(egi.to_vertices_and_faces()) == pytest.approx(4 / 3)


def test_egi_from_vectors_tetrahedron():
    egi = egi_from_vectors(TETRAHEDRON, [0, 0, 0])

    assert egi.number_of_vertices() == 4
    assert egi.number_of_edges() == 6
    assert egi.number_of_faces() == 4
    assert egi.is_closed()
    assert volume_polyhedron(egi.to_vertices_and_faces()) > 0
    for vkey in egi.vertices():
        assert egi.vertex_attribute(vkey, "type") == "face"
        assert TOL.is_close(distance_point_point(egi.vertex_coordinates(vkey), [0, 0, 0]), 1.0)


def test_egi_from_vectors_fewer_than_four_vectors_raise():
    with pytest.raises(Exception, match="Four or more vectors"):
        egi_from_vectors({0: [1, 0, 0], 1: [0, 1, 0], 2: [-1, -1, -1]}, [0, 0, 0])


@pytest.mark.parametrize(
    "extra",
    [
        [3, 0, 0],  # the same direction as AXES[0], another length
        [cos(0.03), sin(0.03), 0],  # 0.03 rad from AXES[0]: 1 - dot = 4.5e-4 < TOL.approximation
    ],
)
def test_egi_from_vectors_coincident_vectors_raise(extra):
    vectors = dict(AXES)
    vectors[6] = extra

    with pytest.raises(Exception, match="Coincident vectors"):
        egi_from_vectors(vectors, [0, 0, 0])


def test_egi_from_vectors_tol():
    vectors = dict(AXES)
    vectors[6] = [cos(0.03), sin(0.03), 0]

    egi = egi_from_vectors(vectors, [0, 0, 0], tol=1e-4)

    assert sorted(egi.vertices()) == sorted(vectors)
    assert egi.is_closed()
    assert egi.number_of_vertices() - egi.number_of_edges() + egi.number_of_faces() == 2
    assert volume_polyhedron(egi.to_vertices_and_faces()) > 0


@pytest.mark.parametrize("n", [3, 5, 6, 8])
def test_egi_from_vectors_right_prism_has_no_zero_vertices(n):
    vectors = {i: [cos(2 * pi * i / n), sin(2 * pi * i / n), 0] for i in range(n)}
    vectors[n] = [0, 0, 1]
    vectors[n + 1] = [0, 0, -1]

    egi = egi_from_vectors(vectors, [0, 0, 0])

    assert [egi.vertex_attribute(vkey, "type") for vkey in egi.vertices()] == ["face"] * (n + 2)
    assert egi.number_of_vertices() == n + 2
    assert egi.number_of_edges() == 3 * n
    assert egi.number_of_faces() == 2 * n
    assert egi.is_closed()
    assert volume_polyhedron(egi.to_vertices_and_faces()) > 0

    cell = cell_from_egi(egi)

    assert cell.number_of_faces() == n + 2
    assert cell.number_of_vertices() == 2 * n
    assert cell.is_closed()


@pytest.mark.parametrize("random_vectors", [5, 8, 11, 14], indirect=True)
def test_egi_from_vectors_zero_vertices(random_vectors):
    origin = [0.5, 0.5, -1.0]

    egi = egi_from_vectors(random_vectors, origin)

    zeros = [vkey for vkey in egi.vertices() if egi.vertex_attribute(vkey, "type") == "zero"]
    assert zeros
    assert sorted(set(egi.vertices()) - set(zeros)) == sorted(random_vectors)
    for vkey in zeros:
        # a crossing of two (or more) great-circle arcs, on the unit sphere
        assert vkey > max(random_vectors)
        assert egi.vertex_attribute(vkey, "magnitude") == 0
        assert egi.vertex_degree(vkey) >= 4
        xyz = egi.vertex_coordinates(vkey)
        assert TOL.is_close(distance_point_point(xyz, origin), 1.0)
        assert TOL.is_allclose(egi.vertex_attribute(vkey, "normal"), normalize_vector(subtract_vectors(xyz, origin)))

    cell = cell_from_egi(egi)

    assert sorted(fkey for fkey in cell.faces() if cell.face_attribute(fkey, "type") == "zero") == sorted(zeros)


@pytest.mark.parametrize("random_vectors", range(4, 17), indirect=True)
def test_egi_from_vectors_random_valid_sets(random_vectors):
    origin = [1.0, -2.0, 0.5]

    egi = egi_from_vectors(random_vectors, origin)

    assert egi.number_of_vertices() - egi.number_of_edges() + egi.number_of_faces() == 2
    assert egi.is_closed()
    assert volume_polyhedron(egi.to_vertices_and_faces()) > 0
    for fkey in egi.faces():
        assert dot_vectors(egi.face_normal(fkey), subtract_vectors(egi.face_centroid(fkey), origin)) > 0
    for vkey in random_vectors:
        assert TOL.is_allclose(egi.vertex_coordinates(vkey), add_vectors(origin, random_vectors[vkey]))
    points = [egi.vertex_coordinates(vkey) for vkey in egi.vertices()]
    for point in points:
        assert TOL.is_close(distance_point_point(point, origin), 1.0)
    # every crossing is one vertex: no two vertices are the same point on the unit sphere
    assert min(distance_point_point(a, b) for a, b in combinations(points, 2)) > TOL.angular


@pytest.mark.parametrize("index", range(len(FORCE_GROUPS)))
def test_egi_from_vectors_equilibrated_force_vectors(index):
    origin = FORCE_GROUPS[index]["origin"]
    vectors = dict(enumerate(FORCE_GROUPS[index]["vectors"]))

    egi = egi_from_vectors(vectors, origin)

    assert egi.number_of_vertices() - egi.number_of_edges() + egi.number_of_faces() == 2
    assert egi.is_closed()
    assert volume_polyhedron(egi.to_vertices_and_faces()) > 0
    for fkey in egi.faces():
        assert dot_vectors(egi.face_normal(fkey), subtract_vectors(egi.face_centroid(fkey), origin)) > 0


@pytest.mark.parametrize(
    "vectors",
    [
        # a cone: all vectors in the open half-space z > 0
        {0: [1, 0, 1], 1: [0, 1, 1], 2: [-1, 0, 1], 3: [0, -1, 1]},
        {0: [1, 0, 0.2], 1: [0, 1, 0.5], 2: [-1, 0.3, 0.1], 3: [0.2, -1, 0.4], 4: [0.3, 0.3, 1], 5: [-0.4, -0.5, 0.6]},
        # the closed half-space z >= 0: the origin is on the boundary of the hull
        {0: [1, 0, 0], 1: [-1, 0, 0], 2: [0, 1, 0], 3: [0, -1, 0], 4: [0, 0, 1]},
        pytest.param(
            {i: [cos(2 * pi * i / 5), sin(2 * pi * i / 5), 0] for i in range(5)},
            marks=pytest.mark.xfail(
                strict=True,
                reason=(
                    "egi.py:218-222: vectors in the xy-plane give a flat EGI of 2 faces whose dot(normal, centroid - origin) is exactly 0; "
                    "neither face counts as outward, so the cycles are flipped and no error is raised"
                ),
            ),
            id="planar",
        ),
    ],
)
def test_egi_from_vectors_not_spanning_3d_raise(vectors):
    with pytest.raises(Exception, match="not consistently oriented"):
        egi_from_vectors(vectors, [0, 0, 0])


# ==============================================================================
# cell_from_egi
# ==============================================================================


def test_cell_from_egi_octahedron_is_cube():
    egi = egi_from_vectors(AXES, [0, 0, 0])

    cell = cell_from_egi(egi)

    assert type(cell) is Cell
    assert cell.number_of_vertices() == 8
    assert cell.number_of_edges() == 12
    assert cell.number_of_faces() == 6
    assert cell.is_closed()
    # vertex keys are egi face keys, face keys are egi vertex keys
    assert sorted(cell.vertices()) == sorted(egi.faces())
    assert sorted(cell.faces()) == sorted(egi.vertices())
    for vkey in cell.vertices():
        assert TOL.is_allclose(cell.vertex_coordinates(vkey), egi.face_center(vkey))
    for fkey in cell.faces():
        assert cell.face_attribute(fkey, "type") == "face"
        assert TOL.is_allclose(cell.face_normal(fkey), AXES[fkey])
        for vkey in cell.face_vertices(fkey):
            assert distance_point_plane(cell.vertex_coordinates(vkey), (cell.face_centroid(fkey), cell.face_normal(fkey))) < TOL.absolute
    for edge in cell.edges():
        assert TOL.is_close(cell.edge_length(edge), 2 / 3)
    # the vertices are the centroids of the octahedron faces, (+-1/3, +-1/3, +-1/3)
    assert volume_polyhedron(cell.to_vertices_and_faces()) == pytest.approx(8 / 27)


def test_cell_from_egi_tetrahedron():
    egi = egi_from_vectors(TETRAHEDRON, [0, 0, 0])

    cell = cell_from_egi(egi)

    assert cell.number_of_vertices() == 4
    assert cell.number_of_edges() == 6
    assert cell.number_of_faces() == 4
    assert cell.is_closed()
    assert volume_polyhedron(cell.to_vertices_and_faces()) > 0
    for fkey in cell.faces():
        assert TOL.is_allclose(cell.face_normal(fkey), normalize_vector(TETRAHEDRON[fkey]))
