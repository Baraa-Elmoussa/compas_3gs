from collections.abc import Callable
from collections.abc import Sequence
from typing import Any

from compas.geometry import add_vectors
from compas.geometry import dot_vectors
from compas.geometry import scale_vector
from compas_3gs.algorithms.other import golden_section_search
from compas_3gs.diagrams import Cell
from compas_3gs.operations import cell_relocate_face
from compas_3gs.operations import cell_split_indet_face_vertices

__all__ = ["cell_arearise_face"]


def cell_arearise_face(cell: Cell, fkey: int, target_area: float, tol: float = 1e-3, callback: Callable | None = None, callback_args: Any = None) -> None:
    """Move a face of a cell along its normal until it reaches a target area.

    The indeterminate vertices of the face are split first. The distance to move the face is found
    with a golden section search on the difference between the face area and the target area.

    Parameters
    ----------
    cell
        A cell.
    fkey
        Identifier of the face.
    target_area
        The target area of the face.
    tol
        Convergence tolerance of the search, an area: the difference between the face area and the target area.
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

    """
    if callback:
        if not callable(callback):
            raise Exception("Callback is not callable.")

    cell_split_indet_face_vertices(cell, fkey)

    area = cell.face_area(fkey)
    normal = cell.face_normal(fkey)
    center = cell.face_centroid(fkey)

    # --------------------------------------------------------------------------
    #   get move direction
    # --------------------------------------------------------------------------
    sign = 1  # if it needs to get larger...
    if area - target_area > 0:
        sign = -1  # if it needs to be smaller...

    if area * target_area < 0:  # if it needs to flip...
        sign = -1

    move_dir = _get_move_direction(cell, fkey) * sign

    # --------------------------------------------------------------------------
    #   evaluation function
    # --------------------------------------------------------------------------

    def evaluation(x):
        xyz = add_vectors(center, scale_vector(normal, x * move_dir))
        new_area = _evaluate_new_face_area(cell, fkey, xyz, normal)

        # callback / conduit ---------------------------------------------------
        if callback:
            callback(cell, callback_args)

        return abs(new_area - target_area)

    # --------------------------------------------------------------------------
    #   golden section search
    # --------------------------------------------------------------------------
    a = 0
    b = abs(area - target_area)

    z = golden_section_search(evaluation, a, b, tol=tol)

    # --------------------------------------------------------------------------
    #   update cell
    # --------------------------------------------------------------------------
    xyz = add_vectors(center, scale_vector(normal, z * move_dir))
    cell_relocate_face(cell, fkey, xyz, normal)


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   helpers
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


def _get_move_direction(cell: Cell, fkey: int) -> int:

    normal = cell.face_normal(fkey)
    area = cell.face_area(fkey)
    center = cell.face_center(fkey)

    new_center = add_vectors(center, normal)

    new_area = _evaluate_new_face_area(cell, fkey, new_center)

    if new_area > area:
        move_dir = 1

    if new_area < area:
        move_dir = -1

    cell_relocate_face(cell, fkey, center, normal)

    if new_area == area:
        raise Exception("Pulling this face will not change its area!")

    return move_dir


def _evaluate_new_face_area(cell: Cell, fkey: int, xyz: Sequence[float], init_normal: Sequence[float] | None = None) -> float:

    if init_normal is None:
        init_normal = cell.face_normal(fkey)

    cell_relocate_face(cell, fkey, xyz, init_normal)

    new_area = cell.face_area(fkey)
    new_normal = cell.face_normal(fkey)

    if dot_vectors(init_normal, new_normal) < 0:
        new_area *= -1

    return new_area


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
