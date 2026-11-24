from math import cos
from math import radians
from math import sin
from math import sqrt

import pytest

from compas.geometry import add_vectors
from compas.geometry import centroid_points
from compas.geometry import distance_point_point
from compas.geometry import normalize_vector
from compas.geometry import scale_vector
from compas.geometry import subtract_vectors
from compas.tolerance import TOL
from compas_3gs.algorithms import EquilibriumResult
from compas_3gs.algorithms import truss_equilibrium

APEX = [0.0, 0.0, 3.0]
FEET = [[2.0 * cos(radians(a)), 2.0 * sin(radians(a)), 0.0] for a in (0, 120, 240)]


def test_tripod_forces_and_reactions():
    vertices = [APEX] + FEET
    edges = [(0, 1), (0, 2), (0, 3)]
    load = [1.5, 0.0, -10.0]
    result = truss_equilibrium(vertices, edges, {0: load}, {1: "xyz", 2: "xyz", 3: "xyz"})

    assert isinstance(result, EquilibriumResult)
    assert result.determinate
    assert TOL.is_zero(result.residual)
    assert all(force < 0 for force in result.forces)  # the legs carry the load in compression
    # the apex is in equilibrium: the member forces push it away from the feet
    total = list(load)
    for (i, j), force in zip(edges, result.forces):
        total = add_vectors(total, scale_vector(normalize_vector(subtract_vectors(vertices[j], vertices[i])), force))
    assert TOL.is_allclose(total, [0, 0, 0])
    # the reactions balance the load
    reaction = [sum(r[k] for r in result.reactions.values()) for k in range(3)]
    assert TOL.is_allclose(add_vectors(reaction, load), [0, 0, 0])


def test_space_truss_on_pinned_slider_and_roller():
    xyz = [[2.0 * cos(radians(a)), 2.0 * sin(radians(a)), 0.0] for a in (90, 210, 330)]
    side = distance_point_point(xyz[0], xyz[1])
    xyz.append(add_vectors(centroid_points(xyz), [0.0, 0.0, side * sqrt(2.0 / 3.0)]))
    for k in range(4, 9):
        xyz.append(subtract_vectors(scale_vector(add_vectors(add_vectors(xyz[k - 1], xyz[k - 2]), xyz[k - 3]), 2.0 / 3.0), xyz[k - 4]))
    edges = [(0, 1), (1, 2), (2, 0)] + [(k, k - i) for k in range(3, 9) for i in (1, 2, 3)]
    loads = {k: [0.5, 0.0, -2.0] for k in range(3, 9)}
    supports = {0: "xyz", 1: "xz", 2: "z"}

    result = truss_equilibrium(xyz, edges, loads, supports)

    assert result.determinate
    assert TOL.is_zero(result.residual)
    assert result.reactions[1][1] == 0.0  # the slider is free along y
    assert result.reactions[2][0] == 0.0 and result.reactions[2][1] == 0.0  # the roller carries z only
    total = [sum(load[k] for load in loads.values()) + sum(r[k] for r in result.reactions.values()) for k in range(3)]
    assert TOL.is_allclose(total, [0, 0, 0])


def test_indeterminate_raises():
    corners = [[1.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [-1.0, -1.0, 0.0], [1.0, -1.0, 0.0]]
    with pytest.raises(ValueError, match="indeterminate"):
        truss_equilibrium([[0.0, 0.0, 2.0]] + corners, [(0, 1), (0, 2), (0, 3), (0, 4)], {0: [0, 0, -10]}, {i: "xyz" for i in range(1, 5)})


def test_mechanism_raises_for_a_load_it_cannot_carry():
    vertices = [[0.0, 0.0, 2.0], [1.0, 1.0, 0.0], [-1.0, -1.0, 0.0]]
    with pytest.raises(ValueError, match="mechanism"):
        truss_equilibrium(vertices, [(0, 1), (0, 2)], {0: [1.0, -1.0, -10.0]}, {1: "xyz", 2: "xyz"})


def test_mechanism_carries_a_load_in_its_plane():
    vertices = [[0.0, 0.0, 2.0], [1.0, 1.0, 0.0], [-1.0, -1.0, 0.0]]
    result = truss_equilibrium(vertices, [(0, 1), (0, 2)], {0: [0.0, 0.0, -10.0]}, {1: "xyz", 2: "xyz"})
    assert not result.determinate
    assert TOL.is_close(result.forces[0], result.forces[1])
    assert TOL.is_close(result.forces[0], -5.0 * distance_point_point(vertices[0], vertices[1]) / 2.0)


def test_tension_is_positive():
    result = truss_equilibrium([[0.0, 0.0, 0.0], [0.0, 0.0, -1.0]], [(0, 1)], {1: [0.0, 0.0, -3.0]}, {0: "xyz"})
    assert TOL.is_close(result.forces[0], 3.0)
    assert TOL.is_allclose(result.reactions[0], [0.0, 0.0, 3.0])
