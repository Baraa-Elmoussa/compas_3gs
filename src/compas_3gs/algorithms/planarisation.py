from collections.abc import Callable
from collections.abc import Sequence
from typing import Any

from compas.datastructures import VolMesh
from compas.geometry import centroid_points
from compas.geometry import distance_point_point
from compas.geometry import dot_vectors
from compas.geometry import project_point_plane
from compas_3gs.datastructures import Mesh3gs
from compas_3gs.diagrams import Cell
from compas_3gs.operations import cell_collapse_short_edge
from compas_3gs.utilities import print_result
from compas_3gs.utilities import scale_polygon

__all__ = ["cell_planarise", "volmesh_planarise"]


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   mesh planarisation
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


def cell_planarise(
    cell: Cell,
    kmax: int = 100,
    target_centers: dict[int, Sequence[float]] | None = None,
    target_normals: dict[int, Sequence[float]] | None = None,
    target_areas: dict[int, float] | None = None,
    fix_vkeys: list[int] | None = None,
    avg_fkeys: list[int] | None = None,
    collapse_edge_length: float = 0.1,
    tolerance_flat: float = 0.001,
    tolerance_area: float = 0.001,
    callback: Callable | None = None,
    callback_args: Any = None,
    print_result_info: bool = False,
) -> None:
    """Planarise the faces of a cell.

    Planarisation of a cell is implemented as a three-step iterative procedure.
    At every iteration, each face is first individually projected to its best-fit plane (unless a target normal is given).
    Then, each face is re-sized to its target area (if given).
    Finally, the new vertex coordinates are computed by taking the centroid of the disconnected corners of the faces.

    Parameters
    ----------
    cell
        A cell.
    kmax
        Number of iterations.
    target_centers
        A dictionary of fkeys and target face centers.
    target_normals
        A dictionary of fkeys and target face normals.
    target_areas
        A dictionary of fkeys and target areas.
    fix_vkeys
        List of vkeys to omit from arearisation.
    avg_fkeys
        List of fkeys to average areas.
    collapse_edge_length
        Minimum length of edge to collapse.
    tolerance_flat
        Convergence tolerance for face flatness, a length: the largest distance of a vertex to its projection plane.
    tolerance_area
        Convergence tolerance for face areas, an area: the largest difference between a face area and its target.
    callback
        A user-defined callback function to be executed after every iteration,
        called as `callback(cell, k, callback_args)`.
    callback_args
        Additional parameters to be passed to the callback.

    Raises
    ------
    Exception
        If a callback is provided, but it is not callable.

    """
    if callback:
        if not callable(callback):
            raise Exception("Callback is not callable.")

    target_centers = target_centers or {}
    target_normals = target_normals or {}
    target_areas = target_areas or {}
    fix_vkeys = fix_vkeys or []
    avg_fkeys = avg_fkeys or []

    # --------------------------------------------------------------------------
    #   1. initialise
    # --------------------------------------------------------------------------
    free_vkeys = list(set(cell.vertices()) - set(fix_vkeys))

    # --------------------------------------------------------------------------
    #   2. loop
    # --------------------------------------------------------------------------
    for k in range(kmax):
        deviation_flat = 0
        deviation_area = 0
        deviation_perp = 0

        new_xyz = {vkey: [] for vkey in cell.vertices()}

        # faces to be averaged -------------------------------------------------
        if avg_fkeys:
            avg_face_area = _avg_face_area(cell, avg_fkeys)

        for fkey in cell.faces():
            # evaluate current face --------------------------------------------
            f_vkeys = cell.face_vertices(fkey)
            f_v_xyz = cell.face_coordinates(fkey)
            f_normal = cell.face_normal(fkey)
            f_area = cell.face_area(fkey)
            f_center = cell.face_centroid(fkey)

            if fkey in target_centers:
                f_center = target_centers[fkey]

            if fkey in target_normals:
                target_normal = target_normals[fkey]

                # perpness deviation
                perpness = 1 - abs(dot_vectors(f_normal, target_normal))
                if perpness > deviation_perp:
                    deviation_perp = perpness

                f_normal = target_normal

            # projection plane -------------------------------------------------
            plane = (f_center, f_normal)

            # ------------------------------------------------------------------
            #   3. planarise
            # ------------------------------------------------------------------
            projected_face = []
            for xyz in f_v_xyz:
                projected_xyz = project_point_plane(xyz, plane)
                projected_face.append(projected_xyz)

                # planarisation deviation
                dist = distance_point_point(xyz, projected_xyz)
                if dist > deviation_flat:
                    deviation_flat = dist

            # ------------------------------------------------------------------
            #   4. arearise
            # ------------------------------------------------------------------
            if fkey in target_areas:
                target_area = target_areas[fkey]

            else:
                target_area = f_area

            if fkey in avg_fkeys:
                target_area = avg_face_area

            # scale factor -----------------------------------------------------
            if target_area != 0:
                scale = (target_area / f_area) ** 0.5

            elif target_area == 0:
                scale = 1 - f_area * 0.1
                # scale = (target_area / f_area) ** 0.5
                # scale = 0.9

            # scale ------------------------------------------------------------
            scaled_face = scale_polygon(projected_face, scale)

            # arearisation deviation
            areaness = abs(f_area - target_area)
            if areaness > deviation_area:
                deviation_area = areaness

            # collect new coordinates
            for i in range(len(f_vkeys)):
                new_xyz[f_vkeys[i]].append(scaled_face[i])

        # ----------------------------------------------------------------------
        #   5. compute new cell vertex coordinates
        # ----------------------------------------------------------------------
        for vkey in free_vkeys:
            final_xyz = centroid_points(new_xyz[vkey])
            cell.vertex_update_xyz(vkey, final_xyz)

        for u, v in cell.edges():
            cell_collapse_short_edge(cell, u, v, min_length=collapse_edge_length)

        # ----------------------------------------------------------------------
        #   7. check convergence
        # ----------------------------------------------------------------------
        if deviation_flat < tolerance_flat and deviation_area < tolerance_area:
            if print_result_info:
                name = "Cell planarisation"
                deviation = deviation_flat

                if target_areas:
                    name = "Cell arearisation"
                    deviation = deviation_area

                print_result(name, k, deviation)

            break

        # callback / conduit ---------------------------------------------------
        if callback:
            callback(cell, k, callback_args)


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   volmesh planarisation
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


def volmesh_planarise(
    volmesh: VolMesh,
    kmax: int = 100,
    target_centers: dict[int, Sequence[float]] | None = None,
    target_normals: dict[int, Sequence[float]] | None = None,
    target_areas: dict[int, float] | None = None,
    fix_vkeys: list[int] | None = None,
    fix_boundary_normals: bool = False,
    fix_all_normals: bool = False,
    min_area: float | None = None,
    max_area: float | None = None,
    tolerance_flat: float = 0.001,
    tolerance_area: float = 0.001,
    tolerance_perp: float = 0.001,
    callback: Callable | None = None,
    callback_args: Any = None,
    print_result_info: bool = False,
) -> None:
    """Planarises the halffaces of a volmesh.

    Planarisation of a volmesh is implemented as a three-step iterative procedure.
    At every iteration, each face is first individually projected to its best-fit plane (unless a target normal is given).
    Then, each face is re-sized to its target area (if given).
    Finally, the new vertex coordinates are computed by taking the centroid of the disconnected corners of the faces.
    Each face is handled once, for both of its halffaces; a target can be keyed by either halfface.

    Parameters
    ----------
    volmesh
        A volmesh object.
    kmax
        Number of iterations.
    target_centers
        A dictionary of halfface keys and target face centers.
    target_normals
        A dictionary of halfface keys and target face normals.
    target_areas
        A dictionary of halfface keys and target areas.
    fix_vkeys
        List of vertices that do not move.
    fix_boundary_normals
        Whether to keep the initial normals of the boundary faces.
    fix_all_normals
        Whether to keep the initial normals of all faces.
    min_area
        Not used.
    max_area
        Not used.
    tolerance_flat
        Convergence tolerance for face flatness, a length: the largest distance of a vertex to its projection plane.
    tolerance_area
        Convergence tolerance for face areas, an area: the largest difference between a face area and its target.
    tolerance_perp
        Convergence tolerance for the face normals: the largest `1 - abs(dot(n, t))` of a unit normal `n`
        and its unit target `t`.
    callback
        A user-defined callback function to be executed after every iteration,
        called as `callback(volmesh, k, callback_args)`.
    callback_args
        Additional parameters to be passed to the callback.
    print_result_info
        If True, print the result of the algorithm.

    Raises
    ------
    Exception
        If a callback is provided, but it is not callable.

    """
    if callback:
        if not callable(callback):
            raise Exception("Callback is not callable.")

    target_centers = target_centers or {}
    target_normals = target_normals or {}
    target_areas = target_areas or {}
    fix_vkeys = fix_vkeys or []

    # --------------------------------------------------------------------------
    #   1. initialise
    # --------------------------------------------------------------------------
    free_vkeys = list(set(volmesh.vertices()) - set(fix_vkeys))
    initial_normals = _get_current_volmesh_normals(volmesh)
    boundary_fkeys = volmesh.halffaces_on_boundaries()

    # --------------------------------------------------------------------------
    #   2. loop
    # --------------------------------------------------------------------------
    for k in range(kmax):
        deviation_flat = 0
        deviation_area = 0
        deviation_perp = 0

        new_xyz = {vkey: [] for vkey in volmesh.vertices()}

        for fkey in volmesh.faces():
            fkey_pair = volmesh.halfface_opposite_halfface(fkey)

            # evaluate current face --------------------------------------------
            f_vkeys = volmesh.halfface_vertices(fkey)
            f_v_xyz = volmesh.halfface_coordinates(fkey)
            f_center = volmesh.halfface_center(fkey)
            f_normal = volmesh.halfface_normal(fkey)
            f_area = volmesh.halfface_area(fkey)

            # override with manual target values -------------------------------
            # (keyed by this halfface or by its opposite; both define the same plane)
            target_key = fkey if fkey in target_centers else fkey_pair
            if target_key in target_centers:
                f_center = target_centers[target_key]

            target_key = fkey if fkey in target_normals else fkey_pair
            if target_key in target_normals:
                target_normal = target_normals[target_key]

                # perpness deviation
                perpness = 1 - abs(dot_vectors(f_normal, target_normal))
                if perpness > deviation_perp:
                    deviation_perp = perpness

                f_normal = target_normal

            if fix_boundary_normals:
                if fkey in boundary_fkeys:
                    f_normal = initial_normals[fkey]["normal"]
            if fix_all_normals:
                f_normal = initial_normals[fkey]["normal"]

            # projection plane -------------------------------------------------
            plane = (f_center, f_normal)

            # ------------------------------------------------------------------
            #   3. planarise
            # ------------------------------------------------------------------
            new_face = []
            for xyz in f_v_xyz:
                projected_xyz = project_point_plane(xyz, plane)
                new_face.append(projected_xyz)

                # planarisation deviation
                flatness = distance_point_point(xyz, projected_xyz)
                if flatness > deviation_flat:
                    deviation_flat = flatness

            # ------------------------------------------------------------------
            #   4. arearise
            # ------------------------------------------------------------------
            if target_areas:
                target_key = fkey if fkey in target_areas else fkey_pair
                if target_key in target_areas:
                    target_area = target_areas[target_key]
                    scale = (target_area / f_area) ** 0.5

                    # scale
                    new_face = scale_polygon(new_face, scale)

                    # arearisation deviation
                    areaness = abs(f_area - target_area)
                    if areaness > deviation_area:
                        deviation_area = areaness

            # collect new coordinates ------------------------------------------
            for i in range(len(f_vkeys)):
                new_xyz[f_vkeys[i]].append(new_face[i])

        # ----------------------------------------------------------------------
        #   5. compute new volmesh vertex coordinates
        # ----------------------------------------------------------------------
        for vkey in free_vkeys:
            final_xyz = centroid_points(new_xyz[vkey])
            volmesh.vertex_update_xyz(vkey, final_xyz)

        # ----------------------------------------------------------------------
        #   6. check convergence
        # ----------------------------------------------------------------------
        if deviation_flat < tolerance_flat and deviation_area < tolerance_area and deviation_perp < tolerance_perp:
            if print_result_info:
                name = "Volmesh planarisation"
                deviation = deviation_flat

                if target_areas:
                    name = "Volmesh arearisation"
                    deviation = deviation_area

                print_result(name, k, deviation)

            break

        # callback / conduit ---------------------------------------------------
        if callback:
            callback(volmesh, k, callback_args)


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   planarisation helpers
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


def _get_current_volmesh_normals(volmesh: VolMesh) -> dict[int, dict]:
    normal_dict = {}
    for hfkey in volmesh.halffaces():
        center = volmesh.halfface_center(hfkey)
        normal = volmesh.halfface_normal(hfkey)
        normal_dict[hfkey] = {"normal": normal, "center": center}
    return normal_dict


def _avg_face_area(mesh: Mesh3gs, fkeys: list[int]) -> float:
    total_area = sum(mesh.face_area(fkey) for fkey in fkeys)
    return total_area / len(fkeys)


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
