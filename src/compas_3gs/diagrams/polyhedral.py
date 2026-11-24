from typing import Any

from compas.geometry import angle_vectors
from compas_3gs.datastructures import Network3gs
from compas_3gs.datastructures import VolMesh3gs

__all__ = ["FormNetwork", "ForceVolMesh"]


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   polyheral form diagram - as network
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


class FormNetwork(Network3gs):
    """A polyhedral form diagram, represented as a network object."""

    __module__ = "compas_3gs.diagrams"

    def __init__(self, *args: Any, name: str = "FormNetwork", **kwargs: Any) -> None:
        super().__init__(*args, name=name, **kwargs)
        self.dual = None

        self.attributes.update(
            {
                "scale": 1,
            }
        )
        self.default_node_attributes.update(
            {
                "is_fixed": False,
                "x_fix": False,
                "y_fix": False,
                "z_fix": False,
            }
        )
        self.default_edge_attributes.update(
            {
                "l_min": 1,
                "l_max": 20,
                "target_vector": None,
                "target_length": None,
                "_a": 0.0,
            }
        )

    def dual_halfface(self, edge: tuple[int, int]) -> int:
        """Find the halfface of the dual force diagram that corresponds to an edge.

        Parameters
        ----------
        edge
            The edge `(u, v)`.

        Returns
        -------
        int
            The halfface of cell `u` that is shared with cell `v` in `self.dual`.

        """
        u, v = edge
        u_halfface, v_halfface = self.dual.cell_pair_halffaces(u, v)
        return u_halfface

    def update_angle_deviations(self) -> None:
        """Compute the angle deviation with the corresponding halfface normal in the ForceVolMesh.

        The angle, in degrees, is stored as the attribute `_a` of each edge and of its halfface.

        """
        for edge in self.edges():
            halfface = self.dual_halfface(edge)
            normal = self.dual.halfface_normal(halfface)
            edge_vector = self.edge_direction(edge)
            a = angle_vectors(edge_vector, normal, deg=True)
            self.edge_attribute(edge, "_a", a)
            self.dual.face_attribute(halfface, "_a", a)


# ******************************************************************************
# ******************************************************************************
# ******************************************************************************
#
#   polyhedral force diagram
#
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************


class ForceVolMesh(VolMesh3gs):
    """A polyhedral force diagram, represented as a volmesh object."""

    __module__ = "compas_3gs.diagrams"

    def __init__(self, *args: Any, name: str = "ForceVolMesh", **kwargs: Any) -> None:
        super().__init__(*args, name=name, **kwargs)
        self.primal = None

        self.attributes.update(
            {
                "scale": 1,
                "convention": -1,
            }
        )
        self.default_vertex_attributes.update(
            {
                "x_fix": False,
                "y_fix": False,
                "z_fix": False,
                "is_fixed": False,
            }
        )
        self.default_edge_attributes.update(
            {
                "l_min": 1,
                "l_max": 20,
                "target_vector": None,
                "target_length": None,
            }
        )
        self.default_face_attributes.update(
            {
                "target_area": None,
                "target_normal": None,
                "_a": 0.0,
            }
        )
        self.default_cell_attributes.update(
            {
                "dir": None,
            }
        )


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
