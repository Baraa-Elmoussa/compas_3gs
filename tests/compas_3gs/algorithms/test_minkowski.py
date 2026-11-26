import random
from math import cos
from math import pi
from math import radians
from math import sin

import pytest

from compas.geometry import angle_vectors
from compas.geometry import centroid_points
from compas.geometry import dot_vectors
from compas.geometry import length_vector
from compas.geometry import normalize_vector
from compas.geometry import subtract_vectors
from compas.geometry import volume_polyhedron
from compas_3gs.algorithms import cell_from_forces
from compas_3gs.algorithms import cell_from_normals_and_areas
from compas_3gs.algorithms import truss_equilibrium
from compas_3gs.diagrams import Cell

AXES = {0: [1, 0, 0], 1: [-1, 0, 0], 2: [0, 1, 0], 3: [0, -1, 0], 4: [0, 0, 1], 5: [0, 0, -1]}
TETRAHEDRON = {0: [1, 1, 1], 1: [1, -1, -1], 2: [-1, 1, -1], 3: [-1, -1, 1]}


@pytest.fixture
def forces_sets():
    # sets of forces in equilibrium: random forces of very different sizes, and the one that closes them
    rng = random.Random(3)
    sets = []
    for n in (5, 8, 12, 16):
        forces = {}
        for key in range(n - 1):
            size = rng.uniform(0.1, 10.0)
            forces[key] = [size * rng.gauss(0, 1) for _ in range(3)]
        forces[n - 1] = [-sum(force[k] for force in forces.values()) for k in range(3)]
        sets.append(forces)
    return sets


# ==============================================================================
# cell_from_normals_and_areas
# ==============================================================================


def test_box():
    cell = cell_from_normals_and_areas(AXES, {0: 6.0, 1: 6.0, 2: 2.0, 3: 2.0, 4: 3.0, 5: 3.0})

    assert isinstance(cell, Cell)
    assert set(cell.faces()) == set(AXES)
    assert cell.number_of_vertices() == 8
    # faces of 6 = b c, 2 = a c, 3 = a b: a box of 1 by 3 by 2
    xs, ys, zs = zip(*[cell.vertex_coordinates(vertex) for vertex in cell.vertices()])
    assert max(xs) - min(xs) == pytest.approx(1.0)
    assert max(ys) - min(ys) == pytest.approx(3.0)
    assert max(zs) - min(zs) == pytest.approx(2.0)
    assert volume_polyhedron(cell.to_vertices_and_faces()) == pytest.approx(6.0)


def test_face_keys_origin_and_orientation():
    normals = {10 + key: normal for key, normal in TETRAHEDRON.items()}
    cell = cell_from_normals_and_areas(normals, {key: 2.0 for key in normals}, origin=[1.0, 2.0, 3.0])

    assert set(cell.faces()) == set(normals)
    assert centroid_points([cell.vertex_coordinates(vertex) for vertex in cell.vertices()]) == pytest.approx([1.0, 2.0, 3.0])
    assert cell.is_closed()
    assert volume_polyhedron(cell.to_vertices_and_faces()) > 0
    for face, normal in normals.items():
        assert angle_vectors(cell.face_normal(face), normal, deg=True) < 1e-6
        assert cell.face_area(face) == pytest.approx(2.0, rel=1e-6)


@pytest.mark.parametrize(
    "normals, areas, message",
    [
        ({0: [1, 0, 0], 1: [-1, 0, 0], 2: [0, 1, 0]}, {0: 1.0, 1: 1.0, 2: 1.0}, "four faces"),
        (AXES, {0: 1.0, 1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0, 5: 0.0}, "positive"),
        ({0: [1, 0, 0], 1: [2, 0, 0], 2: [-1, 0, 0], 3: [0, 1, 0], 4: [0, -1, 0]}, {key: 1.0 for key in range(5)}, "same normal"),
        ({key: [cos(2 * pi * key / 5), sin(2 * pi * key / 5), 0] for key in range(5)}, {key: 1.0 for key in range(5)}, "one plane"),
        (AXES, {0: 1.0, 1: 2.0, 2: 1.0, 3: 1.0, 4: 1.0, 5: 1.0}, "do not close"),
    ],
    ids=["three-faces", "zero-area", "same-normal", "planar", "open"],
)
def test_raises(normals, areas, message):
    with pytest.raises(ValueError, match=message):
        cell_from_normals_and_areas(normals, areas)


# ==============================================================================
# cell_from_forces
# ==============================================================================


def test_cell_from_forces_reaches_every_force(forces_sets):
    for forces in forces_sets:
        cell = cell_from_forces(forces, [0.0, 0.0, 0.0])

        assert set(cell.faces()) == set(forces)
        assert cell.is_closed()
        assert volume_polyhedron(cell.to_vertices_and_faces()) > 0
        for face, force in forces.items():
            assert cell.face_area(face) == pytest.approx(length_vector(force), rel=1e-6)
            assert angle_vectors(cell.face_normal(face), force, deg=True) < 1e-5
        # convex: every vertex is on the inner side of every face plane
        for face in cell.faces():
            normal, point = cell.face_normal(face), cell.face_centroid(face)
            for vertex in cell.vertices():
                assert dot_vectors(subtract_vectors(cell.vertex_coordinates(vertex), point), normal) < 1e-9


def test_cell_from_forces_with_forces_in_one_plane_of_very_different_sizes():
    # four forces in the xy plane and two small ones out of it: the faces that touch follow from the sizes
    forces = {0: [10.0, 0.0, 0.0], 1: [-6.0, 0.0, 0.0], 2: [0.0, 8.0, 0.0], 3: [-2.0, -3.0, 0.0], 4: [0.0, 0.0, 0.5], 5: [-2.0, -5.0, -0.5]}
    cell = cell_from_forces(forces, [0.0, 0.0, 0.0])

    assert cell.is_closed()
    for face, force in forces.items():
        assert cell.face_area(face) == pytest.approx(length_vector(force), rel=1e-6)
        assert angle_vectors(cell.face_normal(face), force, deg=True) < 1e-5


def test_cell_from_forces_with_areas_a_hundred_times_apart():
    # the forces at a foot of a lattice tower: a leg and a face diagonal of about 15000 kN, the reaction, and two
    # members of about 150 and 500 kN; the first step of the search leaves the polytope when it is not kept inside
    forces = {
        0: [-145.07497139367595, 11.437070331436004, -18.78144500701343],
        1: [-7159.36192061388, 7159.36192061388, -11756.783713022149],
        2: [-38.834357397982004, 492.5993392834764, -63.772043601775],
        3: [0.0, 0.0, -786.3004187220452],
        4: [7343.271249405537, -7663.398330228798, 12625.637620353],
    }
    cell = cell_from_forces(forces, [0.0, 0.0, 0.0])

    assert cell.is_closed()
    for face, force in forces.items():
        assert cell.face_area(face) == pytest.approx(length_vector(force), rel=1e-6)
        assert angle_vectors(cell.face_normal(face), force, deg=True) < 1e-5


def test_cell_from_forces_tripod_apex():
    apex = [0.0, 0.0, 3.0]
    feet = [[2.0 * cos(radians(a)), 2.0 * sin(radians(a)), 0.0] for a in (0, 120, 240)]
    load = [1.5, 0.0, -10.0]
    result = truss_equilibrium([apex] + feet, [(0, 1), (0, 2), (0, 3)], {0: load}, {1: "xyz", 2: "xyz", 3: "xyz"})
    forces = {i: [force * c for c in normalize_vector(subtract_vectors(foot, apex))] for i, (foot, force) in enumerate(zip(feet, result.forces))}
    forces[3] = load

    cell = cell_from_forces(forces, apex)

    assert cell.number_of_faces() == 4
    assert centroid_points([cell.vertex_coordinates(vertex) for vertex in cell.vertices()]) == pytest.approx(apex)
    for face, force in forces.items():
        assert cell.face_area(face) == pytest.approx(length_vector(force), rel=1e-6)


def test_cell_from_forces_leaves_out_zero_forces():
    forces = dict(TETRAHEDRON)
    forces[7] = [0.0, 0.0, 0.0]
    cell = cell_from_forces(forces, [0, 0, 0])
    assert set(cell.faces()) == set(TETRAHEDRON)


@pytest.mark.parametrize(
    "forces, message",
    [
        ({0: [1, 0, 0], 1: [-1, 0, 0], 2: [0, 1, 0], 3: [0, 0, 0]}, "four or more"),
        ({0: [1, 0, 0], 1: [2, 0, 0], 2: [0, 1, 0], 3: [0, -1, 0], 4: [-3, 0, 0]}, "same direction"),
        ({i: [cos(2 * pi * i / 5), sin(2 * pi * i / 5), 0] for i in range(5)}, "one plane"),
        ({0: [1, 0, 0], 1: [-1, 0, 0], 2: [0, 1, 0], 3: [0, -1, 0], 4: [0, 0, 1], 5: [0, 0, -2]}, "do not close"),
    ],
    ids=["three-forces", "same-direction", "planar", "not-in-equilibrium"],
)
def test_cell_from_forces_raises(forces, message):
    with pytest.raises(ValueError, match=message):
        cell_from_forces(forces, [0, 0, 0])
