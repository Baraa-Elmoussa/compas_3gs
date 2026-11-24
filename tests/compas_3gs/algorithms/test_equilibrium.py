from math import cos
from math import radians
from math import sin
from math import sqrt

import numpy as np
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


PYRAMID = [[0.2, -0.1, 2.0], [1.0, 1.0, 0.0], [-1.0, 1.5, 0.0], [-1.0, -1.0, 0.0], [1.2, -1.0, 0.0]]
LEGS = [(0, 1), (0, 2), (0, 3), (0, 4)]
PINNED = {i: "xyz" for i in range(1, 5)}
APEX_LOAD = {0: [1.0, 0.5, -10.0]}


def test_indeterminate_truss_has_one_self_stress_state():
    result = truss_equilibrium(PYRAMID, LEGS, APEX_LOAD, PINNED)

    assert not result.determinate
    assert len(result.self_stresses) == 1
    member_forces, reactions = result.self_stresses[0]
    assert max(abs(force) for force in member_forces) == pytest.approx(1.0)
    # a self-stress state is in equilibrium without any load
    total = [0.0, 0.0, 0.0]
    for (i, j), force in zip(LEGS, member_forces):
        total = add_vectors(total, scale_vector(normalize_vector(subtract_vectors(PYRAMID[j], PYRAMID[i])), force))
    assert TOL.is_allclose(total, [0, 0, 0])
    for (i, j), force in zip(LEGS, member_forces):
        foot = add_vectors(reactions[j], scale_vector(normalize_vector(subtract_vectors(PYRAMID[i], PYRAMID[j])), force))
        assert TOL.is_allclose(foot, [0, 0, 0])


def test_indeterminate_elastic_state_is_the_finite_element_state():
    ea = [1.0, 2.0, 3.0, 4.0]
    result = truss_equilibrium(PYRAMID, LEGS, APEX_LOAD, PINNED, stiffness=ea)

    # the apex is the only free node: its stiffness matrix is the sum of EA / L * d d^T over the legs
    K = np.zeros((3, 3))
    for (i, j), stiffness in zip(LEGS, ea):
        d = np.subtract(PYRAMID[j], PYRAMID[i])
        K += stiffness / np.linalg.norm(d) ** 3 * np.outer(d, d)
    u = np.linalg.solve(K, APEX_LOAD[0])
    fem = []
    for (i, j), stiffness in zip(LEGS, ea):
        d = np.subtract(PYRAMID[j], PYRAMID[i])
        length = np.linalg.norm(d)
        fem.append(stiffness / length * np.dot(-d / length, u))  # the foot j does not move, the apex i moves by u

    assert TOL.is_allclose(result.forces, fem)
    assert not TOL.is_allclose(result.forces, truss_equilibrium(PYRAMID, LEGS, APEX_LOAD, PINNED).forces)


def test_indeterminate_state_from_a_prescribed_force_or_reaction():
    by_force = truss_equilibrium(PYRAMID, LEGS, APEX_LOAD, PINNED, forces={2: -1.5})
    assert TOL.is_close(by_force.forces[2], -1.5)
    assert TOL.is_zero(by_force.residual)

    by_reaction = truss_equilibrium(PYRAMID, LEGS, APEX_LOAD, PINNED, reactions={(1, "z"): 3.0})
    assert TOL.is_close(by_reaction.reactions[1][2], 3.0)
    assert TOL.is_zero(by_reaction.residual)

    # one prescribed value fixes the one free parameter, whatever the stiffness
    stiff = truss_equilibrium(PYRAMID, LEGS, APEX_LOAD, PINNED, forces={2: -1.5}, stiffness=[1.0, 10.0, 0.1, 5.0])
    assert TOL.is_allclose(stiff.forces, by_force.forces)


def test_prescribed_values_that_cannot_hold_raise():
    with pytest.raises(ValueError, match="cannot all hold"):
        truss_equilibrium([APEX] + FEET, [(0, 1), (0, 2), (0, 3)], {0: [0, 0, -10]}, {1: "xyz", 2: "xyz", 3: "xyz"}, forces={0: 5.0})
    with pytest.raises(ValueError, match="cannot all hold"):
        truss_equilibrium(PYRAMID, LEGS, APEX_LOAD, PINNED, forces={0: 0.0, 1: 0.0})
    with pytest.raises(ValueError, match="no reaction"):
        truss_equilibrium(PYRAMID, LEGS, APEX_LOAD, PINNED, reactions={(0, "z"): 1.0})
    with pytest.raises(ValueError, match="positive"):
        truss_equilibrium(PYRAMID, LEGS, APEX_LOAD, PINNED, stiffness=[1.0, 0.0, 1.0, 1.0])


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
