from collections.abc import Sequence
from typing import Any

from compas.datastructures import VolMesh
from compas.geometry import centroid_points
from compas.geometry import normalize_vector

__all__ = ["VolMesh3gs"]


class VolMesh3gs(VolMesh):
    """Inherits and extends the compas VolMesh class, such that it is more suitable for 3D graphic statics applications.

    Primarily used for polyhedral form and force diagrams.

    """

    def __init__(self, *args: Any, name: str = "VolMesh3gs", **kwargs: Any) -> None:
        super().__init__(*args, name=name, **kwargs)
        self.primal = None
        self.attributes.update(
            {
                "scale": 1.0,
            }
        )
        self.default_face_attributes.update(
            {
                "_a": 0.0,
            }
        )

    # --------------------------------------------------------------------------
    #   vertices
    # --------------------------------------------------------------------------

    def vertex_normal(self, vertex: int) -> list[float] | None:
        """Return the normal vector at a boundary vertex as the normalised average of the
        area-weighted normals of its boundary halffaces.

        Parameters
        ----------
        vertex
            The identifier of the vertex.

        Returns
        -------
        list[float] | None
            The components of the normal vector, or `None` if the vertex is not on the boundary.

        """
        if not self.is_vertex_on_boundary(vertex):
            return

        halffaces = []
        for halfface in self.vertex_halffaces(vertex):
            if self.is_halfface_on_boundary(halfface):
                halffaces.append(halfface)

        vectors = [self.face_normal(halfface, False) for halfface in halffaces if halfface is not None]
        return normalize_vector(centroid_points(vectors))

    def vertex_update_xyz(self, vertex: int, new_xyz: Sequence[float], constrained: bool = True) -> None:
        """Update the coordinates of a vertex, respecting its per-axis fixity.

        Parameters
        ----------
        vertex
            The identifier of the vertex.
        new_xyz
            The new XYZ coordinates.
        constrained
            If True, an axis is only updated if the vertex attribute `x_fix`, `y_fix` or `z_fix`
            is not True (a missing attribute counts as False).
            If False, all three coordinates are updated.

        """
        for axis, value in zip("xyz", new_xyz):
            if constrained and self.vertex_attribute(vertex, axis + "_fix"):
                continue
            self.vertex_attribute(vertex, axis, value)

    # --------------------------------------------------------------------------
    # cell
    # --------------------------------------------------------------------------

    def cell_pair_halffaces(self, cell_1: int, cell_2: int) -> tuple[int, int]:
        """Find the pair of halffaces through which two cells touch.

        Parameters
        ----------
        cell_1
            Identifier of the cell 1.
        cell_2
            Identifier of the cell 2.

        Returns
        -------
        tuple[int, int]
            The halfface of cell 1 and the opposite halfface, which belongs to cell 2.

        Raises
        ------
        KeyError
            If the two cells do not share a face.

        """
        for halfface in self.cell_faces(cell_1):
            if self.halfface_opposite_cell(halfface) == cell_2:
                return halfface, self.halfface_opposite_halfface(halfface)

        raise KeyError((cell_1, cell_2))


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
