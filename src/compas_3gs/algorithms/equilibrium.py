from collections.abc import Sequence
from typing import NamedTuple

import numpy as np

from compas.geometry import normalize_vector
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
        True if the truss is stable, so that it carries any load; False for a mechanism that carries this load only.

    """

    forces: list[float]
    reactions: dict[int, list[float]]
    residual: float
    determinate: bool


def truss_equilibrium(
    vertices: Sequence[Sequence[float]],
    edges: Sequence[tuple[int, int]],
    loads: dict[int, Sequence[float]],
    supports: dict[int, str],
    tol: float | None = None,
) -> EquilibriumResult:
    """Compute the axial forces and support reactions of a pin-jointed truss from equilibrium.

    Every vertex gives three equilibrium equations. The unknowns are one axial force per edge and
    one reaction component per restrained direction of every support. The forces follow from statics
    alone only if the equations have a single solution: no more unknowns than independent equations.

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
    tol
        The largest out-of-balance force, relative to the largest load, for the truss to carry the load.
        Default is `TOL.absolute`.

    Returns
    -------
    EquilibriumResult

    Raises
    ------
    ValueError
        If the truss is statically indeterminate (its forces depend on the member stiffness),
        or if it is a mechanism that cannot carry the load.

    Examples
    --------
    >>> vertices = [[0, 0, 1], [-1, 0, 0], [1, 0, 0], [0, 1, 0]]
    >>> result = truss_equilibrium(vertices, [(0, 1), (0, 2), (0, 3)], {0: [0, 0, -2]}, {1: "xyz", 2: "xyz", 3: "xyz"})
    >>> [round(force, 3) + 0.0 for force in result.forces]
    [-1.414, -1.414, 0.0]
    >>> [round(component, 3) + 0.0 for component in result.reactions[1]]
    [1.0, 0.0, 1.0]

    """
    n = len(vertices)
    columns = []
    for i, j in edges:
        column = np.zeros(3 * n)
        direction = normalize_vector(subtract_vectors(vertices[j], vertices[i]))  # a tension pulls i towards j
        column[3 * i : 3 * i + 3] += direction
        column[3 * j : 3 * j + 3] -= direction
        columns.append(column)
    restrained = []
    for vertex, dofs in supports.items():
        for axis in dofs:
            column = np.zeros(3 * n)
            column[3 * vertex + AXES.index(axis)] = 1.0
            columns.append(column)
            restrained.append((vertex, AXES.index(axis)))
    A = np.array(columns).reshape(-1, 3 * n).T
    p = np.zeros(3 * n)
    for vertex, load in loads.items():
        p[3 * vertex : 3 * vertex + 3] += load

    unknowns = A.shape[1]
    rank = int(np.linalg.matrix_rank(A)) if unknowns else 0
    if rank < unknowns:
        raise ValueError(
            "The truss is statically indeterminate (degree {}): {} unknowns ({} members and {} reaction components) "
            "but {} independent equilibrium equations, so the forces depend on the member stiffness.".format(unknowns - rank, unknowns, len(edges), len(restrained), rank)
        )
    x = np.linalg.lstsq(A, -p, rcond=None)[0] if unknowns else np.zeros(0)
    residual = float(np.max(np.abs(A @ x + p))) if n else 0.0
    scale = max(1.0, float(np.max(np.abs(p)))) if n else 1.0
    if not TOL.is_zero(residual / scale, tol):
        raise ValueError("The truss is a mechanism that cannot carry the load: {} equations of rank {}, out of balance by {:.3g}.".format(3 * n, rank, residual))

    reactions = {vertex: [0.0, 0.0, 0.0] for vertex in supports}
    for (vertex, axis), value in zip(restrained, x[len(edges) :]):
        reactions[vertex][axis] = float(value) + 0.0
    return EquilibriumResult([float(force) for force in x[: len(edges)]], reactions, residual, rank == 3 * n)
