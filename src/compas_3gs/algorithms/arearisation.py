from collections.abc import Callable
from collections.abc import Sequence
from typing import Any

from compas.geometry import add_vectors
from compas.geometry import dot_vectors
from compas.geometry import scale_vector
from compas.tolerance import TOL
from compas_3gs.algorithms.other import golden_section_search
from compas_3gs.diagrams import Cell
from compas_3gs.operations import cell_relocate_face
from compas_3gs.operations import cell_split_indet_face_vertices

__all__ = ["cell_arearise_face"]


def cell_arearise_face(
    cell: Cell,
    fkey: int,
    target_area: float,
    tol: float = 1e-3,
    zmax: float | None = None,
    callback: Callable | None = None,
    callback_args: Any = None,
) -> None:
    """Move a face of a cell along its normal until it reaches a target area.

    The indeterminate vertices of the face are split first; the other faces keep their planes. The area of the face
    is a parabola in its position z along its normal. A golden section search finds the z where the face has the
    target area, on the range from 0 up to `zmax` or up to the vertex of the parabola, whichever comes first: of two
    positions with the target area, the nearer one.

    Parameters
    ----------
    cell
        A cell.
    fkey
        Identifier of the face.
    target_area
        The target area of the face. A negative target area turns the face over: its normal ends opposite to its
        initial normal.
    tol
        Convergence tolerance of the search, a distance along the normal: the search stops when its range is shorter.
    zmax
        The end of the search range, a distance along the normal. Default is the difference between the area of the
        face and the target area.
    callback
        A user-defined callback function to be executed after every evaluation,
        called as `callback(cell, callback_args)`.
    callback_args
        Additional parameters to be passed to the callback.

    Raises
    ------
    Exception
        If a callback is provided, but it is not callable,
        or if moving the face along its normal does not change its area.
    ValueError
        If the face does not reach the target area within the search range.

    """
    if callback:
        if not callable(callback):
            raise Exception("Callback is not callable.")

    cell_split_indet_face_vertices(cell, fkey)

    area = cell.face_area(fkey)
    normal = cell.face_normal(fkey)
    center = cell.face_centroid(fkey)
    face = {vertex: cell.vertex_coordinates(vertex) for vertex in cell.face_vertices(fkey)}

    # --------------------------------------------------------------------------
    #   the area as a parabola of the position along the normal
    # --------------------------------------------------------------------------
    h = max(cell.face_area(other) for other in cell.faces()) ** 0.5
    area_plus = _face_area_at(cell, fkey, face, center, normal, h)
    area_minus = _face_area_at(cell, fkey, face, center, normal, -h)
    _face_area_at(cell, fkey, face, center, normal, 0.0)

    if TOL.is_close(area_plus, area) and TOL.is_close(area_minus, area):
        raise Exception("Pulling this face will not change its area!")

    # A(z) = area + slope * z + curvature * z**2
    slope = (area_plus - area_minus) / (2 * h)
    curvature = (area_plus + area_minus - 2 * area) / (2 * h * h)

    # --------------------------------------------------------------------------
    #   move direction: the one in which the area goes towards the target
    # --------------------------------------------------------------------------
    grow = 1 if target_area > area else -1
    if TOL.is_close(area_plus, area_minus):
        move_dir = 1
    else:
        move_dir = grow if slope > 0 else -grow
    slope *= move_dir

    # --------------------------------------------------------------------------
    #   search range: up to zmax, and no further than the vertex of the parabola
    # --------------------------------------------------------------------------
    b = abs(area - target_area) if zmax is None else zmax
    if curvature and 0 < -slope / (2 * curvature) < b:
        b = -slope / (2 * curvature)

    area_b = area + slope * b + curvature * b * b
    if (area_b - target_area) * (area - target_area) > 0:
        raise ValueError("Face {} goes from area {:.4g} to {:.4g} over a move of {:.4g} along its normal, and does not reach {:.4g}.".format(fkey, area, area_b, b, target_area))

    # --------------------------------------------------------------------------
    #   golden section search
    # --------------------------------------------------------------------------

    def evaluation(z):
        new_area = _face_area_at(cell, fkey, face, center, normal, z * move_dir)

        # callback / conduit ---------------------------------------------------
        if callback:
            callback(cell, callback_args)

        return abs(new_area - target_area)

    z = golden_section_search(evaluation, 0, b, tol=tol)

    # --------------------------------------------------------------------------
    #   update cell
    # --------------------------------------------------------------------------
    _face_area_at(cell, fkey, face, center, normal, z * move_dir)


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   helpers
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


def _face_area_at(cell: Cell, fkey: int, face: dict[int, Sequence[float]], center: Sequence[float], normal: Sequence[float], z: float) -> float:
    # move the face from its initial vertices `face`, so that every vertex slides along its initial edge line,
    # to the plane through center + z * normal, and return its area, negative if the face turned over
    for vertex, xyz in face.items():
        cell.vertex_attributes(vertex, "xyz", xyz)
    cell_relocate_face(cell, fkey, add_vectors(center, scale_vector(normal, z)), normal)

    area = cell.face_area(fkey)
    if dot_vectors(normal, cell.face_normal(fkey)) < 0:
        area *= -1
    return area


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   Main
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


if __name__ == "__main__":
    pass
