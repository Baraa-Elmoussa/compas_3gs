from collections.abc import Sequence
from typing import NamedTuple

import numpy as np

from compas.geometry import length_vector
from compas.geometry import scale_vector
from compas.geometry import subtract_vectors
from compas.tolerance import TOL

__all__ = ["EquilibriumResult", "truss_equilibrium"]

AXES = "xyz"


class EquilibriumResult(NamedTuple):
    """The member forces and support reactions of a truss in equilibrium.

    Attributes
    ----------
    forces
        The axial force of each edge, in the order of the input edges, tension positive.
    reactions
        The reaction vector of each support vertex.
    residual
        The largest out-of-balance force component at any vertex.
    determinate
        True if the equilibrium equations have one solution for any load: the truss has no self-stress state and no mechanism.
    self_stresses
        The self-stress states of the truss, each as member forces and reactions in the format of `forces` and `reactions`,
        scaled so that the largest member force is 1. Adding any combination of them to a state in equilibrium gives
        another state in equilibrium with the same loads. Empty for a statically determinate truss.

    """

    forces: list[float]
    reactions: dict[int, list[float]]
    residual: float
    determinate: bool
    self_stresses: list[tuple[list[float], dict[int, list[float]]]]


def truss_equilibrium(
    vertices: Sequence[Sequence[float]],
    edges: Sequence[tuple[int, int]],
    loads: dict[int, Sequence[float]],
    supports: dict[int, str],
    stiffness: Sequence[float] | None = None,
    forces: dict[int, float] | None = None,
    reactions: dict[tuple[int, str], float] | None = None,
    tol: float | None = None,
) -> EquilibriumResult:
    """Compute the axial forces and support reactions of a pin-jointed truss in equilibrium.

    Every vertex gives three equilibrium equations. The unknowns are one axial force per edge and one reaction component
    per restrained direction of every support. If the equations have one solution, the truss is statically determinate
    and statics alone gives the forces. If they have more, the truss is statically indeterminate: the states in equilibrium
    are one of them plus any combination of the self-stress states, one free parameter per state.

    Of those states, the function returns the one with the least complementary energy, the sum of `N**2 * L / (2 * EA)`
    over the members, among the states that have the prescribed `forces` and `reactions`. With the stiffness `EA` of the
    real members and rigid supports, and nothing prescribed, this is the state of the linear elastic truss, the one a
    linear finite element model gives. Prescribing as many independent forces or reactions as there are self-stress states
    fixes the state whatever the stiffness.

    Parameters
    ----------
    vertices
        The vertex coordinates.
    edges
        The edges, as pairs of vertex indices.
    loads
        The force applied to each loaded vertex.
    supports
        The restrained directions of each support vertex, as a string of axes,
        for example `"xyz"` for a pinned support or `"z"` for a roller.
    stiffness
        The axial stiffness `EA` of each edge. Default is the same for all edges.
    forces
        Prescribed axial forces, by edge index.
    reactions
        Prescribed reaction components, by support vertex and axis, for example `{(4, "z"): 10.0}`.
    tol
        The largest out-of-balance force, relative to the largest load, for the truss to carry the load.
        Default is `TOL.absolute`.

    Returns
    -------
    EquilibriumResult

    Raises
    ------
    ValueError
        If the truss is a mechanism that cannot carry the load, if a stiffness is not positive,
        if a prescribed reaction is not one of the restrained directions, or if the prescribed forces and reactions
        cannot all hold in equilibrium.

    Examples
    --------
    >>> vertices = [[0, 0, 1], [-1, 0, 0], [1, 0, 0], [0, 1, 0]]
    >>> result = truss_equilibrium(vertices, [(0, 1), (0, 2), (0, 3)], {0: [0, 0, -2]}, {1: "xyz", 2: "xyz", 3: "xyz"})
    >>> [round(force, 3) + 0.0 for force in result.forces]
    [-1.414, -1.414, 0.0]
    >>> [round(component, 3) + 0.0 for component in result.reactions[1]]
    [1.0, 0.0, 1.0]

    Four legs carry a load at their apex: one leg too many for statics, so one self-stress state.

    >>> vertices = [[0, 0, 1], [1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]
    >>> edges = [(0, 1), (0, 2), (0, 3), (0, 4)]
    >>> supports = {1: "xyz", 2: "xyz", 3: "xyz", 4: "xyz"}
    >>> result = truss_equilibrium(vertices, edges, {0: [0, 0, -4]}, supports)
    >>> result.determinate, len(result.self_stresses)
    (False, 1)
    >>> [round(force, 3) for force in result.forces]
    [-1.414, -1.414, -1.414, -1.414]
    >>> result = truss_equilibrium(vertices, edges, {0: [0, 0, -4]}, supports, forces={0: 0.0})
    >>> [round(force, 3) + 0.0 for force in result.forces]
    [0.0, -2.828, 0.0, -2.828]

    """
    n = len(vertices)
    m = len(edges)
    columns = []
    lengths = []
    for i, j in edges:
        vector = subtract_vectors(vertices[j], vertices[i])
        length = length_vector(vector)
        direction = scale_vector(vector, 1 / length)  # a tension pulls i towards j
        column = np.zeros(3 * n)
        column[3 * i : 3 * i + 3] += direction
        column[3 * j : 3 * j + 3] -= direction
        columns.append(column)
        lengths.append(length)
    restrained = []
    for vertex, dofs in supports.items():
        for axis in dofs:
            column = np.zeros(3 * n)
            column[3 * vertex + AXES.index(axis)] = 1.0
            columns.append(column)
            restrained.append((vertex, axis))
    A = np.array(columns).reshape(-1, 3 * n).T
    p = np.zeros(3 * n)
    for vertex, load in loads.items():
        p[3 * vertex : 3 * vertex + 3] += load

    # the states in equilibrium: one solution x0 plus any combination of the self-stress states, the null space of A
    unknowns = A.shape[1]
    rank = 0
    states = np.zeros((unknowns, 0))
    x = np.zeros(unknowns)
    if unknowns and n:
        _, s, vt = np.linalg.svd(A)
        rank = int(np.sum(s > s.max() * max(A.shape) * np.finfo(float).eps)) if s.max() > 0 else 0
        states = vt[rank:].T
        x = np.linalg.lstsq(A, -p, rcond=None)[0]
    residual = float(np.max(np.abs(A @ x + p))) if n else 0.0
    scale = max(1.0, float(np.max(np.abs(p)))) if n else 1.0
    if not TOL.is_zero(residual / scale, tol):
        raise ValueError("The truss is a mechanism that cannot carry the load: {} equations of rank {}, out of balance by {:.3g}.".format(3 * n, rank, residual))

    ea = np.ones(m) if stiffness is None else np.array(stiffness, dtype=float)
    if np.any(ea <= 0):
        raise ValueError("The stiffness of every member must be positive.")
    rows = list((forces or {}).keys())
    values = list((forces or {}).values())
    for (vertex, axis), value in (reactions or {}).items():
        if (vertex, axis) not in restrained:
            raise ValueError("Vertex {} has no reaction along {}.".format(vertex, axis))
        rows.append(m + restrained.index((vertex, axis)))
        values.append(value)

    k = states.shape[1]
    if k:
        # least complementary energy, with the prescribed values as constraints (a KKT system)
        weights = np.array(lengths) / ea
        H = states[:m].T @ (weights[:, None] * states[:m])
        g = states[:m].T @ (weights * x[:m])
        C = states[rows]
        kkt = np.block([[H, C.T], [C, np.zeros((len(rows), len(rows)))]])
        rhs = np.concatenate([-g, np.array(values, dtype=float) - x[rows]])
        x = x + states @ np.linalg.lstsq(kkt, rhs, rcond=None)[0][:k]
    if rows and not TOL.is_allclose(x[rows], values):
        raise ValueError("The prescribed forces and reactions cannot all hold: statics fixes some of them, or they are more than the {} self-stress states allow.".format(k))

    # the state and each self-stress state as member forces and reactions by vertex
    results = []
    for state in [x] + [state / state[:m][np.argmax(np.abs(state[:m]))] for state in states.T]:
        components = {vertex: [0.0, 0.0, 0.0] for vertex in supports}
        for (vertex, axis), value in zip(restrained, state[m:]):
            components[vertex][AXES.index(axis)] = float(value) + 0.0
        results.append(([float(force) + 0.0 for force in state[:m]], components))
    return EquilibriumResult(results[0][0], results[0][1], residual, k == 0 and rank == 3 * n, results[1:])
