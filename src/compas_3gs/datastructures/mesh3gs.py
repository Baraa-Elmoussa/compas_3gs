from collections.abc import Sequence
from typing import Any

from compas.datastructures import Mesh
from compas.geometry import Vector
from compas.geometry import cross_vectors
from compas.geometry import length_vector
from compas.geometry import normal_polygon
from compas.geometry import normalize_vector
from compas.geometry import subtract_vectors
from compas.tolerance import TOL

__all__ = ["Mesh3gs"]


class Mesh3gs(Mesh):
    """Inherits and extends the compas Mesh class, such that it is more suitable for 3D graphic statics applications.

    Primarily used for the EGI.

    """

    def __init__(self, *args: Any, name: str = "Mesh3gs", **kwargs: Any) -> None:
        super().__init__(*args, name=name, **kwargs)

    # --------------------------------------------------------------------------
    # helpers - vertices
    # --------------------------------------------------------------------------

    def vertex_update_xyz(self, vkey: int, xyz: Sequence[float], constrained: bool = True) -> None:
        """Update the coordinates of a vertex, respecting its per-axis fixity.

        Parameters
        ----------
        vkey
            The identifier of the vertex.
        xyz
            The new XYZ coordinates.
        constrained
            If True, an axis is only updated if the vertex attribute `x_fix`, `y_fix` or `z_fix`
            is not True (a missing attribute counts as False).
            If False, all three coordinates are updated.

        """
        for axis, value in zip("xyz", xyz):
            if constrained and self.vertex_attribute(vkey, axis + "_fix"):
                continue
            self.vertex_attribute(vkey, axis, value)

    # --------------------------------------------------------------------------
    # helpers - faces
    # --------------------------------------------------------------------------

    def face_area(self, fkey: int) -> float:
        """Compute the vector area of a face.

        Parameters
        ----------
        fkey
            The identifier of the face.

        Returns
        -------
        float
            The length of the area-weighted normal of the face.

        """
        return length_vector(normal_polygon(self.face_coordinates(fkey), unitized=False))

    def face_normal(self, fkey: int, unitized: bool = True) -> Vector:
        """Compute the normal of a face.

        Parameters
        ----------
        fkey
            The identifier of the face.
        unitized
            If True, the normal is unitized.

        Returns
        -------
        Vector
            The normal of the face. For a face with a zero vector area, the unitized cross product
            of its first two edges, which is always a unit vector.

        """
        points = self.face_coordinates(fkey)
        normal = normal_polygon(points, unitized=unitized)
        if TOL.is_zero(length_vector(normal)):
            uv = subtract_vectors(points[1], points[0])
            vw = subtract_vectors(points[2], points[1])
            normal = normalize_vector(cross_vectors(uv, vw))
        return Vector(*normal)


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
