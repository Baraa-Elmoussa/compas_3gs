from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from compas.datastructures import Mesh

from compas.geometry import subtract_vectors
from compas.geometry import normalize_vector
from compas.geometry import normal_polygon
from compas.geometry import length_vector
from compas.geometry import cross_vectors


__all__ = ['Mesh3gs']


class Mesh3gs(Mesh):
    """Inherits and extends the compas Mesh class, such that it is more suitable for 3D graphic statics applications.

    Primarily used for the EGI.

    """

    def __init__(self, *args, name="Mesh3gs", **kwargs):
        super().__init__(*args, name=name, **kwargs)

    # --------------------------------------------------------------------------
    # helpers - vertices
    # --------------------------------------------------------------------------

    def vertex_update_xyz(self, vkey, xyz, constrained=True):
        """Update the coordinates of a vertex, respecting its per-axis fixity.

        Parameters
        ----------
        vkey : int
            The identifier of the vertex.
        xyz : list of float
            The new XYZ coordinates.
        constrained : bool, optional
            If True, an axis is only updated if the vertex attribute ``x_fix``, ``y_fix`` or ``z_fix``
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

    def face_area(self, fkey):
        """Compute the vector area of a face.

        Parameters
        ----------
        fkey : int
            The identifier of the face.

        Returns
        -------
        float
            The length of the area-weighted normal of the face.


        """
        return length_vector(normal_polygon(self.face_coordinates(fkey), unitized=False))

    def face_normal(self, fkey, unitized=True):
        """Compute the normal of a face.

        Parameters
        ----------
        fkey : int
            The identifier of the face.
        unitized : bool, optional
            If True, the normal is unitized.

        Returns
        -------
        list of float
            The normal of the face.

        Notes
        -----
        Equal to :meth:`compas.datastructures.Mesh.face_normal`, except for faces with a zero vector area
        (for example a bowtie). There the normal falls back to the cross product of the first two edges,
        because ``cell_planarise`` uses the face normal as the projection plane and a zero normal
        makes that plane degenerate.

        """
        points = self.face_coordinates(fkey)
        normal = normal_polygon(points, unitized=unitized)
        if length_vector(normal) == 0:
            uv = subtract_vectors(points[1], points[0])
            vw = subtract_vectors(points[2], points[1])
            normal = normalize_vector(cross_vectors(uv, vw))
        return list(normal)


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   Main
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


if __name__ == '__main__':
    pass
