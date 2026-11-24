from collections.abc import Sequence
from typing import Any

from compas.datastructures import Graph

__all__ = ["Network3gs"]


class Network3gs(Graph):
    """Inherits and extends the compas Graph class, such that it is more suitable for 3D graphic statics applications.

    Primarily used for polyhedral (and possibly non-polyhedral) form diagrams.

    """

    def __init__(self, *args: Any, name: str = "Network3gs", **kwargs: Any) -> None:
        super().__init__(*args, name=name, **kwargs)
        self.dual = None
        self.default_edge_attributes.update(
            {
                "lmin": 0.0,
                "lmax": 1e7,
                "_a": 0.0,
            }
        )

    # --------------------------------------------------------------------------
    # helpers - vertices
    # --------------------------------------------------------------------------

    def vertex_update_xyz(self, node: int, new_xyz: Sequence[float], constrained: bool = True) -> None:
        """Update the coordinates of a node, respecting its per-axis fixity.

        Parameters
        ----------
        node
            The identifier of the node.
        new_xyz
            The new XYZ coordinates.
        constrained
            If True, an axis is only updated if the node attribute `x_fix`, `y_fix` or `z_fix`
            is not True (a missing attribute counts as False).
            If False, all three coordinates are updated.

        """
        for axis, value in zip("xyz", new_xyz):
            if constrained and self.node_attribute(node, axis + "_fix"):
                continue
            self.node_attribute(node, axis, value)

    # --------------------------------------------------------------------------
    # helpers - edges
    # --------------------------------------------------------------------------

    def edge_avg_length(self) -> float:
        """Compute the average length of the edges.

        Returns
        -------
        float
            The average edge length.

        """
        return sum(self.edge_length(edge) for edge in self.edges()) / self.number_of_edges()

    # --------------------------------------------------------------------------
    # drawing
    # --------------------------------------------------------------------------

    # def draw(self, **kwattr):
    #     artist = NetworkArtist(self)
    #     artist.draw_edges(**kwattr)
    #     artist.draw_nodes(**kwattr)

    # def clear(self, **kwattr):
    #     artist = NetworkArtist(self)
    #     artist.clear_by_name()
    #     artist.clear_layer()

    # def draw_nodes(self, **kwattr):
    #     artist = NetworkArtist(self)
    #     artist.draw_nodes(**kwattr)

    # def draw_edges(self, **kwattr):
    #     artist = NetworkArtist(self)
    #     artist.draw_edges(**kwattr)

    # def draw_vertexlabels(self, **kwattr):
    #     artist = NetworkArtist(self)
    #     artist.draw_vertexlabels(**kwattr)

    # def draw_edgelabels(self, **kwattr):
    #     artist = NetworkArtist(self)
    #     artist.draw_edgelabels(**kwattr)


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
