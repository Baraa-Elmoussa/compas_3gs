from collections.abc import Sequence

import numpy as np
from scipy.optimize import linprog
from scipy.optimize import minimize
from scipy.spatial import ConvexHull
from scipy.spatial import HalfspaceIntersection

from compas.geometry import cross_vectors
from compas.geometry import dot_vectors
from compas.geometry import length_vector
from compas.geometry import normalize_vector
from compas.tolerance import TOL
from compas_3gs.diagrams import Cell

__all__ = ["cell_from_normals_and_areas", "cell_from_forces"]


def _polytope(normals: np.ndarray, offsets: np.ndarray) -> tuple[np.ndarray, np.ndarray, list[list[int]], float] | None:
    """The polytope `normals @ x <= offsets`: its vertices, the area of each face, the vertices of each face, its volume.

    Returns None if the polytope is empty or flat.

    """
    n = len(normals)
    # a point inside: the centre of the largest ball in the polytope, by linear programming
    ball = linprog(np.array([0.0, 0.0, 0.0, -1.0]), A_ub=np.hstack([normals, np.ones((n, 1))]), b_ub=offsets, bounds=[(None, None)] * 3 + [(0.0, None)])
    if not ball.success or ball.x[3] <= TOL.absolute * max(1.0, np.abs(offsets).max()):
        return None
    points = HalfspaceIntersection(np.hstack([normals, -offsets[:, None]]), ball.x[:3]).intersections
    # several planes through one corner give the corner several times
    size = max(1.0, np.abs(points).max())
    welded = []
    for point in points:
        if all(np.linalg.norm(point - other) > TOL.absolute * size for other in welded):
            welded.append(point)
    points = np.array(welded)
    hull = ConvexHull(points)
    areas = np.zeros(n)
    for simplex, equation in zip(hull.simplices, hull.equations):
        a, b, c = points[simplex]
        areas[int(np.argmax(normals @ equation[:3]))] += 0.5 * np.linalg.norm(np.cross(b - a, c - a))
    distance = np.abs(points @ normals.T - offsets[None, :])  # of every vertex to every face plane
    faces = [list(np.flatnonzero(distance[:, i] <= TOL.absolute * size)) for i in range(n)]
    return points, areas, faces, hull.volume


def cell_from_normals_and_areas(
    normals: dict[int, Sequence[float]],
    areas: dict[int, float],
    origin: Sequence[float] | None = None,
    tol: float | None = None,
) -> Cell:
    """Construct the convex cell with the given outward face normals and face areas.

    By Minkowski's theorem, unit normals that span space and positive areas with `sum(area * normal) = 0`
    are the faces of exactly one convex polyhedron, up to a translation. Which faces meet, and along which
    edges, follows from the areas as much as from the normals.

    The cell is found as Minkowski proved it exists: the distances `h` of the face planes from the centre
    that give the largest volume for a fixed `sum(area * h)` give face areas proportional to the target
    areas. A quasi-Newton optimisation finds them (it minimises `sum(area * h) - log(volume) / 3`, a
    convex function of `h`), the cell is scaled to the target areas, and Newton steps on the face areas,
    with the derivative of each face area with respect to the plane distances, bring every area to its
    target.

    Parameters
    ----------
    normals
        The outward normal of each face, by face key; they do not need to be unit vectors.
    areas
        The area of each face, by face key, all positive.
    origin
        The centroid of the vertices of the cell. Default is the world origin.
    tol
        The largest relative difference between a face area and its target, and between `sum(area * normal)`
        and zero relative to the sum of the areas. Default is `TOL.relative`.

    Returns
    -------
    Cell
        The cell, its face keys those of `normals`.

    Raises
    ------
    ValueError
        If there are fewer than four faces, an area is not positive, two normals point the same way, the normals
        lie in one plane, or the faces do not close (`sum(area * normal)` is not zero).

    Examples
    --------
    >>> normals = {0: [1, 0, 0], 1: [-1, 0, 0], 2: [0, 1, 0], 3: [0, -1, 0], 4: [0, 0, 1], 5: [0, 0, -1]}
    >>> areas = {0: 6.0, 1: 6.0, 2: 2.0, 3: 2.0, 4: 3.0, 5: 3.0}
    >>> cell = cell_from_normals_and_areas(normals, areas)
    >>> round(cell.face_area(0), 6), round(cell.face_area(4), 6)
    (6.0, 3.0)
    >>> cell.number_of_vertices()
    8

    """
    tol = TOL.relative if tol is None else tol
    keys = list(normals)
    if len(keys) < 4:
        raise ValueError("A cell needs four faces or more, not {}.".format(len(keys)))
    target = np.array([float(areas[key]) for key in keys])
    if np.any(target <= 0):
        raise ValueError("Every face area must be positive.")
    units = [normalize_vector(normals[key]) for key in keys]
    N = np.array(units)
    for i, a in enumerate(units):
        for b in units[i + 1 :]:
            if dot_vectors(a, b) > 1 - TOL.approximation:
                raise ValueError("Two faces have the same normal.")
    if np.linalg.svd(N, compute_uv=False)[-1] < TOL.angular:
        raise ValueError("The face normals lie in one plane.")
    closure = np.linalg.norm(target @ N) / target.sum()
    if closure > tol:
        raise ValueError("The faces do not close: sum(area * normal) is {:.3g} of the sum of the areas.".format(closure))

    # 1. the shape: minimise sum(w * h) - log(volume) / 3, convex in h, with the areas w scaled to sum to 1;
    #    its gradient w - areas / (3 volume) is zero where the face areas are proportional to w
    weights = target / target.sum()

    def objective(offsets: np.ndarray) -> tuple[float, np.ndarray]:
        polytope = _polytope(N, offsets)
        if polytope is None:
            return np.inf, np.zeros_like(offsets)
        _, faces_area, _, volume = polytope
        return weights @ offsets - np.log(volume) / 3, weights - faces_area / (3 * volume)

    offsets = minimize(objective, np.ones(len(keys)), jac=True, method="L-BFGS-B").x

    # 2. the size, then Newton steps on the areas: dA_i/dh_j = l_ij / sin(t_ij), dA_i/dh_i = -sum_j l_ij cos(t_ij) / sin(t_ij)
    points, faces_area, faces, _ = _polytope(N, offsets)
    offsets = offsets * np.sqrt(target.sum() / faces_area.sum())
    for _ in range(20):
        points, faces_area, faces, _ = _polytope(N, offsets)
        error = faces_area - target
        if np.max(np.abs(error) / target) <= tol:
            break
        J = np.zeros((len(keys), len(keys)))
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                shared = sorted(set(faces[i]) & set(faces[j]))
                if len(shared) < 2:
                    continue
                length = max(np.linalg.norm(points[p] - points[q]) for p in shared for q in shared)
                cos, sin = dot_vectors(units[i], units[j]), length_vector(cross_vectors(units[i], units[j]))
                J[i, j] = J[j, i] = length / sin
                J[i, i] -= length * cos / sin
                J[j, j] -= length * cos / sin
        offsets = offsets + np.linalg.lstsq(J, -error, rcond=None)[0]  # the translations are the null space of J
    if np.max(np.abs(faces_area - target) / target) > tol:
        raise ValueError("No cell reaches the areas: largest relative area error {:.3g}.".format(np.max(np.abs(faces_area - target) / target)))

    # 3. the cell: each face's vertices in turn around its normal
    shift = (np.zeros(3) if origin is None else np.array(origin, dtype=float)) - points.mean(axis=0)
    cell = Cell()
    for index, point in enumerate(points + shift):
        cell.add_vertex(key=index, x=float(point[0]), y=float(point[1]), z=float(point[2]))
    for i, key in enumerate(keys):
        corners = points[faces[i]]
        centre = corners.mean(axis=0)
        u = normalize_vector(cross_vectors(units[i], [1.0, 0.0, 0.0] if abs(units[i][0]) < 0.9 else [0.0, 1.0, 0.0]))
        v = cross_vectors(units[i], u)
        turn = np.arctan2((corners - centre) @ np.array(v), (corners - centre) @ np.array(u))
        cell.add_face([int(faces[i][k]) for k in np.argsort(turn)], fkey=key)
    return cell


def cell_from_forces(forces: dict[int, Sequence[float]], origin: Sequence[float], tol: float | None = None) -> Cell:
    """Construct the force polyhedron of a node from the forces that act on it.

    Each force gives one face of the cell: its normal is the direction of the force and its area is the magnitude
    of the force. Forces in equilibrium close the faces into one convex cell (`cell_from_normals_and_areas`).

    Parameters
    ----------
    forces
        The force vectors on the node, by integer key; the key becomes the face key of the cell.
        Forces of zero length are left out.
    origin
        The location of the node: the centroid of the vertices of the cell.
    tol
        The largest relative difference between a face area and its force, and the largest out-of-balance
        force relative to the sum of the force magnitudes. Default is `TOL.relative`.

    Returns
    -------
    Cell

    Raises
    ------
    ValueError
        If fewer than four forces are not zero, if two forces have the same direction, if all forces lie in one
        plane, or if the forces are not in equilibrium.

    Examples
    --------
    >>> forces = {0: [1, 0, 0], 1: [-1, 0, 0], 2: [0, 1, 0], 3: [0, -1, 0], 4: [0, 0, 2], 5: [0, 0, -2]}
    >>> cell = cell_from_forces(forces, [0, 0, 0])
    >>> cell.number_of_faces()
    6
    >>> round(cell.face_area(4), 6)
    2.0

    """
    vectors = {key: list(force) for key, force in forces.items() if not TOL.is_zero(length_vector(force))}
    if len(vectors) < 4:
        raise ValueError("A force polyhedron needs four or more forces, not {}.".format(len(vectors)))
    directions = [normalize_vector(vector) for vector in vectors.values()]
    for i, a in enumerate(directions):
        for b in directions[i + 1 :]:
            if dot_vectors(a, b) > 1 - TOL.approximation:
                raise ValueError("Two forces have the same direction.")
    if np.linalg.svd(np.array(directions), compute_uv=False)[-1] < TOL.angular:
        raise ValueError("The forces lie in one plane.")
    return cell_from_normals_and_areas(vectors, {key: length_vector(vector) for key, vector in vectors.items()}, origin=origin, tol=tol)
