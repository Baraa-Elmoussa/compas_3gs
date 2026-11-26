from typing import TYPE_CHECKING

from compas.geometry import add_vectors
from compas.geometry import bestfit_plane
from compas.geometry import centroid_points
from compas.geometry import centroid_points_weighted
from compas.geometry import cross_vectors
from compas.geometry import distance_point_point
from compas.geometry import length_vector
from compas.geometry import normalize_vector
from compas.geometry import project_point_plane
from compas.geometry import scale_vector
from compas.geometry import subtract_vectors

if TYPE_CHECKING:
    from compas.datastructures import VolMesh
    from compas_3gs.diagrams import Cell


__all__ = [
    "resultant_vector",
    "polygon_normal_oriented",
    "polygon_area_oriented",
    "polygon_area_footprint",
    "scale_polygon",
    "cell_face_flatness",
    "cell_face_areaness",
    "volmesh_face_flatness",
    "volmesh_face_areaness",
    "datastructure_centroid",
]


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   vectors
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


def resultant_vector(vectors, locations):
    """Computes the resultant vector of a group of vectors.

    Parameters
    ----------
    vectors: dictionary
        vector_key-XYZ (vector components) pairs

    locations: dictionary
        vector_key-XYZ (vector locations) pairs

    Returns
    -------
    vector: tuple
        XYZ components of the resultant vector

    point: tuple
        XYZ coordinates of the location of the resultant vector

    Examples
    --------
    >>> vectors = {0: [0.97, 0.83, 0.83], 1: [0.65, 0.89, 0.34], 2: [0.68, 0.63, 0.97]}
    >>> resultant, location = resultant_vector(vectors, vectors)
    >>> resultant
    [2.3, 2.35, 2.14]
    >>> location
    [0.781, 0.780, 0.736]

    """
    points = []
    weights = []

    for vkey in vectors:
        points.append(locations[vkey])
        weights.append(length_vector(vectors[vkey]))

    resultant_xyz = centroid_points_weighted(points, weights)
    x, y, z = zip(*vectors.values())
    resultant_vector = [sum(x), sum(y), sum(z)]

    return resultant_vector, resultant_xyz


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   polygons
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


def polygon_normal_oriented(polygon, unitized=True):
    """Compute the oriented normal of any closed polygon (can be convex, concave or complex).

    Parameters
    ----------
    polygon : sequence
        The XYZ coordinates of the vertices/corners of the polygon.
        The vertices are assumed to be in order.
        The polygon is assumed to be closed:
        the first and last vertex in the sequence should not be the same.

    Returns
    -------
    list
        The weighted or unitized normal vector of the polygon.

    """
    p = len(polygon)
    assert p > 2, "At least three points required"
    w = centroid_points(polygon)
    normal_sum = (0, 0, 0)

    for i in range(-1, len(polygon) - 1):
        u = polygon[i]
        v = polygon[i + 1]
        uv = subtract_vectors(v, u)
        vw = subtract_vectors(w, v)
        normal = scale_vector(cross_vectors(uv, vw), 0.5)
        normal_sum = add_vectors(normal_sum, normal)

    if not unitized:
        return normal_sum

    return normalize_vector(normal_sum)


def polygon_area_oriented(polygon):
    """Compute the oriented area of a polygon (can be convex or concave).

    Parameters
    ----------
    polygon : sequence
        The XYZ coordinates of the vertices/corners of the polygon.
        The vertices are assumed to be in order.
        The polygon is assumed to be closed:
        the first and last vertex in the sequence should not be the same.

    Returns
    -------
    float
        The area of the polygon.

    """
    return length_vector(polygon_normal_oriented(polygon, unitized=False))


def polygon_area_footprint(polygon):
    """Compute the non-oriented area of a polygon (can be convex or concave).

    Parameters
    ----------
    polygon : list of lists
        A list of polygon point coordinates.

    Returns
    -------
    float
        The non-oriented area of the polygon.

    """
    area = 0
    w = centroid_points(polygon)

    for i in range(-1, len(polygon) - 1):
        u = polygon[i]
        v = polygon[i + 1]
        uv = subtract_vectors(v, u)
        vw = subtract_vectors(w, v)
        normal = scale_vector(cross_vectors(uv, vw), 0.5)
        area += length_vector(normal)

    return area


def scale_polygon(polygon, scale):
    """Scale a polygon.

    Parameters
    ----------
    polygon : list of lists
        A list of polygon point coordinates.

    Returns
    -------
    list of lists
        Reordered polygon point coordinates.

    """
    center = centroid_points(polygon)

    scaled_polygon = []

    for pt in polygon:
        vector = subtract_vectors(pt, center)
        new_pt = add_vectors(center, scale_vector(vector, scale))
        scaled_polygon.append(new_pt)

    return scaled_polygon


def polygon_flatness(polygon):
    """Comput the flatness of a polygon.

    Parameters
    ----------
    polygon : list of lists
        A list of polygon point coordinates.

    Returns
    -------
    float
        The flatness.

    Note
    ----
    compas.geometry.mesh_flatness function currently only works for quadrilateral faces.
    This function uses the distance between each face vertex and its projected point on the best-fit plane of the face as the flatness metric.

    """
    deviation = 0

    plane = bestfit_plane(polygon)

    for pt in polygon:
        pt_proj = project_point_plane(pt, plane)
        dev = distance_point_point(pt, pt_proj)
        if dev > deviation:
            deviation = dev

    return deviation


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   mesh
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


def cell_face_flatness(cell: "Cell") -> dict[int, float]:
    """Compute the flatness of every face of a mesh.

    Parameters
    ----------
    cell
        Cell as a mesh object.

    Returns
    -------
    dict[int, float]
        A dictionary of flatness deviation for every face of the mesh: the largest distance
        of a face vertex to the best-fit plane of the face.

    """

    flatness_dict = {fkey: 0 for fkey in cell.faces()}

    for fkey in cell.faces():
        f_vkeys = cell.face_vertices(fkey)
        f_pts = [cell.vertex_coordinates(vkey) for vkey in f_vkeys]
        dev = polygon_flatness(f_pts)

        flatness_dict[fkey] = dev

    return flatness_dict


def cell_face_areaness(cell: "Cell", target_areas: dict[int, float]) -> dict[int, float]:
    """Compute the areaness of every face of a mesh.

    Parameters
    ----------
    cell
        Cell as a mesh object.
    target_areas
        Dictionary of target areas.

    Returns
    -------
    dict[int, float]
        A dictionary of area deviation for each face of the mesh.

    """
    areaness_dict = {}

    for fkey in target_areas:
        area = cell.face_area(fkey)
        areaness_dict[fkey] = abs(target_areas[fkey] - area)

    return areaness_dict


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   volmesh
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


def volmesh_face_flatness(volmesh: "VolMesh") -> dict[int, float]:
    """Compute the flatness of every face of a volmesh.

    Parameters
    ----------
    volmesh
        A volmesh.

    Returns
    -------
    dict[int, float]
        A dictionary of flatness deviation for every face of the volmesh (one halfface per face):
        the largest distance of a face vertex to the best-fit plane of the face.

    """
    flatness_dict = {fkey: 0 for fkey in volmesh.faces()}

    for fkey in volmesh.faces():
        f_vkeys = volmesh.halfface_vertices(fkey)
        f_pts = [volmesh.vertex_coordinates(vkey) for vkey in f_vkeys]
        dev = polygon_flatness(f_pts)

        flatness_dict[fkey] = dev

    return flatness_dict


def volmesh_face_areaness(volmesh: "VolMesh", target_areas: dict[int, float]) -> dict[int, float]:
    """Compute the areaness of the halffaces of a volmesh that have a target area.

    Parameters
    ----------
    volmesh
        A volmesh.
    target_areas
        Dictionary of halfface keys and target areas.

    Returns
    -------
    dict[int, float]
        A dictionary of area deviation for each halfface in `target_areas`.

    """
    areaness_dict = {}

    for hfkey in target_areas:
        area = volmesh.halfface_area(hfkey)
        areaness_dict[hfkey] = abs(target_areas[hfkey] - area)

    return areaness_dict


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   datastructures
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


def datastructure_centroid(datastructure):
    """Compute the centroid of the vertices of a mesh or volmesh.

    Parameters
    ----------
    datastructure
        A mesh or volmesh object.

    Returns
    -------
    list
        The coordinates of the centroid.

    """
    points = [datastructure.vertex_coordinates(vkey) for vkey in datastructure.vertices()]

    return centroid_points(points)


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
