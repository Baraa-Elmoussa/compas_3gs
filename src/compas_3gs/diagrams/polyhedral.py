from collections.abc import Iterator
from collections.abc import Sequence
from typing import Any

from compas.geometry import Vector
from compas.geometry import angle_vectors
from compas.geometry import distance_point_point
from compas.tolerance import TOL
from compas_3gs.datastructures import Network3gs
from compas_3gs.datastructures import VolMesh3gs

__all__ = ["FormNetwork", "ForceVolMesh"]

AXES = "xyz"
SUPPORT_TYPES = {"pinned": "xyz", "slider": "xz", "roller": "z", "free": ""}


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
    """A polyhedral form diagram, represented as a network object.

    Supports, applied forces and the results of equilibrium are node and edge attributes:
    `support` (the restrained directions, for example `"xyz"`), `px`, `py`, `pz` (the applied force),
    `_rx`, `_ry`, `_rz` (the support reaction) and, per edge, `_f` (the axial force, tension positive).

    Examples
    --------
    >>> form = FormNetwork.from_lines([([0, 0, 1], [-1, 0, 0]), ([0, 0, 1], [1, 0, 0])])
    >>> form.set_support(form.node_at([-1, 0, 0]), "pinned")
    >>> form.set_support(form.node_at([1, 0, 0]), "roller")
    >>> form.add_force(form.node_at([0, 0, 1]), [0, 0, -2])
    >>> [form.node_attribute(node, "support") for node in form.supports()]
    ['xyz', 'z']

    """

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
                "support": "",
                "px": 0.0,
                "py": 0.0,
                "pz": 0.0,
                "_rx": 0.0,
                "_ry": 0.0,
                "_rz": 0.0,
            }
        )
        self.default_edge_attributes.update(
            {
                "l_min": 1,
                "l_max": 20,
                "target_vector": None,
                "target_length": None,
                "_a": 0.0,
                "_f": 0.0,
            }
        )

    def node_at(self, point: Sequence[float], tol: float | None = None) -> int:
        """Find the node at a point.

        Parameters
        ----------
        point
            The point.
        tol
            The largest distance between the point and the node. Default is `TOL.absolute`.

        Returns
        -------
        int

        Raises
        ------
        KeyError
            If there is no node at the point.

        """
        node = min(self.nodes(), key=lambda node: distance_point_point(self.node_coordinates(node), point))
        if not TOL.is_zero(distance_point_point(self.node_coordinates(node), point), tol):
            raise KeyError("There is no node at {}.".format(list(point)))
        return node

    def set_support(self, node: int, dofs: str = "xyz") -> None:
        """Make a node a support, or change its support.

        Parameters
        ----------
        node
            The node.
        dofs
            The restrained directions, as axes (`"xyz"`, `"xz"`, `"z"`, ...) or as a support type:
            `"pinned"` (x, y and z), `"slider"` (x and z), `"roller"` (z) or `"free"` (none).
            The members are pin-jointed and carry no moments, so a fixed support is a pinned support.

        Raises
        ------
        ValueError
            If `dofs` is neither axes nor a support type.

        """
        dofs = SUPPORT_TYPES.get(dofs, dofs)
        if any(axis not in AXES for axis in dofs):
            raise ValueError("The support must be axes from 'xyz' or one of {}, not {!r}.".format(sorted(SUPPORT_TYPES), dofs))
        self.node_attribute(node, "support", "".join(axis for axis in AXES if axis in dofs))

    def add_force(self, node: int, vector: Sequence[float]) -> None:
        """Add a force to a node.

        Parameters
        ----------
        node
            The node.
        vector
            The force vector; it adds to the force already on the node.

        """
        px, py, pz = self.node_attributes(node, ["px", "py", "pz"])
        self.node_attributes(node, ["px", "py", "pz"], [px + vector[0], py + vector[1], pz + vector[2]])

    def supports(self) -> Iterator[int]:
        """Iterate over the nodes with at least one restrained direction.

        Yields
        ------
        int

        """
        for node in self.nodes():
            if self.node_attribute(node, "support"):
                yield node

    def node_load(self, node: int) -> Vector:
        """Return the force applied to a node.

        Parameters
        ----------
        node
            The node.

        Returns
        -------
        Vector

        """
        return Vector(*self.node_attributes(node, ["px", "py", "pz"]))

    def node_reaction(self, node: int) -> Vector:
        """Return the support reaction at a node.

        Parameters
        ----------
        node
            The node.

        Returns
        -------
        Vector

        """
        return Vector(*self.node_attributes(node, ["_rx", "_ry", "_rz"]))

    def edge_force(self, edge: tuple[int, int]) -> float:
        """Return the axial force in a member, tension positive.

        Parameters
        ----------
        edge
            The edge.

        Returns
        -------
        float

        """
        return self.edge_attribute(edge, "_f")

    def node_forces(self, node: int) -> dict[tuple[int, int] | str, Vector]:
        """Return the forces that act on a node: one per member, the applied force and the support reaction.

        A member in tension pulls the node towards its neighbour; a member in compression pushes it away.

        Parameters
        ----------
        node
            The node.

        Returns
        -------
        dict[tuple[int, int] | str, Vector]
            The force of each member, by edge, and, if they are not zero, the applied force (`"load"`) and the reaction (`"reaction"`).

        """
        forces = {}
        for nbr in self.neighbors(node):
            edge = (node, nbr) if self.has_edge((node, nbr)) else (nbr, node)
            direction = Vector.from_start_end(self.node_coordinates(node), self.node_coordinates(nbr)).unitized()
            forces[edge] = direction * self.edge_force(edge)
        load = self.node_load(node)
        if not TOL.is_zero(load.length):
            forces["load"] = load
        reaction = self.node_reaction(node)
        if not TOL.is_zero(reaction.length):
            forces["reaction"] = reaction
        return forces

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
