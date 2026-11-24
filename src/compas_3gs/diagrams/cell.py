from typing import Any

from compas.geometry import centroid_polyhedron
from compas_3gs.datastructures import Mesh3gs

__all__ = ["Cell"]


class Cell(Mesh3gs):
    """A single polyhedral cell, represented as a mesh object."""

    def __init__(self, *args: Any, name: str = "Cell", **kwargs: Any) -> None:
        super().__init__(*args, name=name, **kwargs)

        self.default_vertex_attributes.update(
            {
                "x_fix": False,
                "y_fix": False,
                "z_fix": False,
            }
        )
        self.default_edge_attributes.update(
            {
                "target_vector": None,
                "target_length": None,
            }
        )
        self.default_face_attributes.update(
            {
                "target_area": None,
                "target_normal": None,
            }
        )

    # --------------------------------------------------------------------------
    # misc
    # --------------------------------------------------------------------------

    def cell_center(self) -> list[float]:
        """Compute the volume centroid of the cell.

        Returns
        -------
        list[float]
            The XYZ coordinates of the centroid of the polyhedron.

        """
        return centroid_polyhedron(self.to_vertices_and_faces())


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
