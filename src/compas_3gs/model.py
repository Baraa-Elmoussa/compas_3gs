from collections.abc import Sequence
from typing import Any

from compas.data import Data
from compas.geometry import Vector
from compas.geometry import scale_vector
from compas.tolerance import TOL
from compas_3gs.algorithms import EquilibriumResult
from compas_3gs.algorithms import cell_from_forces
from compas_3gs.algorithms import truss_equilibrium
from compas_3gs.algorithms import volmesh_dual_network
from compas_3gs.algorithms import volmesh_reciprocate
from compas_3gs.diagrams import Cell
from compas_3gs.diagrams import ForceVolMesh
from compas_3gs.diagrams import FormNetwork
from compas_3gs.utilities import get_force_mags

__all__ = ["GSModel"]


class GSModel(Data):
    """A 3D graphic statics model: a form diagram, its equilibrium and its force diagram.

    A model is built from the form diagram (`from_lines`) or from a polyhedral force diagram
    (`from_force_diagram`). Supports and forces go on the nodes of the form diagram; `solve`
    computes the member forces and support reactions from equilibrium, and `derive_force_diagram`
    constructs the force polyhedron of every node. The diagrams stay accessible as `form`, `force`
    and `cells`, for the algorithms and operations of `compas_3gs`.

    Member forces are positive in tension and negative in compression.

    Parameters
    ----------
    form
        The form diagram. Default is an empty `FormNetwork`.
    force
        The polyhedral force diagram of the form diagram, if there is one.
    name
        The name of the model.

    Attributes
    ----------
    form : FormNetwork
        The form diagram, with the supports and forces and, after `solve`, the member forces and reactions.
    force : ForceVolMesh | None
        The polyhedral force diagram.
    result : EquilibriumResult | None
        The result of the last `solve`.
    self_stresses : list[tuple[dict[tuple[int, int], float], dict[int, Vector]]]
        The self-stress states found by the last `solve`, each as member forces by edge and reactions by support node.
        Adding any multiple of one to the forces in the form diagram gives another state in equilibrium with the same loads.
    cells : dict[int, Cell]
        The force polyhedron of each node, after `derive_force_diagram`.
    cell_faces : dict[int, dict[int, tuple[int, int] | str]]
        For each node with a cell, what each face of the cell stands for: an edge of the form diagram, `"load"` or `"reaction"`.
        The other faces of the cell carry no force.
    skipped : dict[int, str]
        The nodes without a force polyhedron, with the reason.

    Examples
    --------
    >>> model = GSModel.from_lines([([0, 0, 1], [-1, 0, 0]), ([0, 0, 1], [1, 0, 0]), ([0, 0, 1], [0, 1, 0])])
    >>> for point in ([-1, 0, 0], [1, 0, 0], [0, 1, 0]):
    ...     model.set_support(model.node_at(point), "pinned")
    >>> model.add_force(model.node_at([0, 0, 1]), [0, 0, -2])
    >>> result = model.solve()
    >>> [round(force, 3) + 0.0 for force in model.forces().values()]
    [-1.414, -1.414, 0.0]

    """

    @property
    def __data__(self) -> dict[str, Any]:
        return {"form": self.form, "force": self.force}

    def __init__(self, form: FormNetwork | None = None, force: ForceVolMesh | None = None, name: str | None = None) -> None:
        super().__init__(name=name)
        self.form = FormNetwork() if form is None else form
        self.force = force
        if force is not None:
            self.form.dual = force
            force.primal = self.form
        self.result = None
        self.self_stresses = []
        self.cells = {}
        self.cell_faces = {}
        self.skipped = {}

    # --------------------------------------------------------------------------
    # constructors
    # --------------------------------------------------------------------------

    @classmethod
    def from_lines(cls, lines: Sequence[tuple[Sequence[float], Sequence[float]]], precision: int | None = None) -> "GSModel":
        """Construct a model from the lines of the members; the lines are joined where their end points meet.

        Parameters
        ----------
        lines
            The start and end point of each member.
        precision
            The precision of the geometric keys that join the end points. Default is `TOL.precision`.

        Returns
        -------
        GSModel

        """
        return cls(FormNetwork.from_lines(lines, precision=precision))

    @classmethod
    def from_force_diagram(cls, force: ForceVolMesh, kmax: int = 500, edge_min: float = 0.1, edge_max: float = 1e6, tolerance: float = 0.001) -> "GSModel":
        """Construct a model from a polyhedral force diagram.

        The form diagram is the dual network of the force diagram, one node per cell and one member per interior face,
        reciprocated so that every member is perpendicular to its face while the force diagram stays as it is.
        The member forces are the face areas, and the boundary faces of each cell are the force applied to its node.

        Parameters
        ----------
        force
            The polyhedral force diagram.
        kmax
            The number of reciprocation iterations.
        edge_min, edge_max
            The shortest and longest member during reciprocation.
        tolerance
            The largest angle, in degrees, between a member and the normal of its face.

        Returns
        -------
        GSModel

        """
        form = volmesh_dual_network(force, cls=FormNetwork)
        volmesh_reciprocate(force, form, kmax=kmax, weight=1.0, edge_min=edge_min, edge_max=edge_max, tolerance=tolerance)
        for edge, value in get_force_mags(force, form).items():
            form.edge_attribute(edge, "_f", value)
        for halfface in force.halffaces_on_boundaries():
            load = scale_vector(force.halfface_normal(halfface), force.halfface_area(halfface) * force.attributes["convention"])
            form.add_force(force.halfface_cell(halfface), load)
        return cls(form, force)

    # --------------------------------------------------------------------------
    # supports and forces
    # --------------------------------------------------------------------------

    def node_at(self, point: Sequence[float], tol: float | None = None) -> int:
        """Find the node of the form diagram at a point.

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
        return self.form.node_at(point, tol=tol)

    def set_support(self, node: int, dofs: str = "xyz") -> None:
        """Make a node of the form diagram a support.

        Parameters
        ----------
        node
            The node.
        dofs
            The restrained directions, as axes (`"xyz"`, `"xz"`, `"z"`, ...) or as a support type:
            `"pinned"` (x, y and z), `"slider"` (x and z), `"roller"` (z) or `"free"` (none).

        """
        self.form.set_support(node, dofs)

    def add_force(self, node: int, vector: Sequence[float]) -> None:
        """Add a force to a node of the form diagram.

        Parameters
        ----------
        node
            The node.
        vector
            The force vector; it adds to the forces already on the node.

        """
        self.form.add_force(node, vector)

    # --------------------------------------------------------------------------
    # analysis
    # --------------------------------------------------------------------------

    def solve(
        self,
        stiffness: float | dict[tuple[int, int], float] | None = None,
        forces: dict[tuple[int, int], float] | None = None,
        reactions: dict[tuple[int, str], float] | None = None,
        tol: float | None = None,
    ) -> EquilibriumResult:
        """Compute member forces and support reactions in equilibrium, and store them in the form diagram.

        A statically determinate truss has one state in equilibrium with its loads. A statically indeterminate truss has
        a family of them, with one free parameter per self-stress state (`EquilibriumResult.self_stresses`). Of those,
        `solve` takes the state with the least complementary energy among the states that have the prescribed `forces`
        and `reactions` (see `truss_equilibrium`): with nothing prescribed, the linear elastic state of the truss with
        the given `stiffness`; with as many prescribed forces or reactions as there are self-stress states, the state
        they fix.

        The member forces go into the edge attribute `_f`, the reactions into the node attributes `_rx`, `_ry` and `_rz`.
        The cells of an earlier `derive_force_diagram` are removed.

        Parameters
        ----------
        stiffness
            The axial stiffness `EA` of the members: one value for all, or one per edge. Default is the same for all.
        forces
            Prescribed axial forces, by edge.
        reactions
            Prescribed reaction components, by support node and axis, for example `{(node, "z"): 10.0}`.
        tol
            The largest out-of-balance force, relative to the largest load; see `truss_equilibrium`.

        Returns
        -------
        EquilibriumResult

        Raises
        ------
        ValueError
            If the truss is a mechanism that cannot carry the load, or if the prescribed forces and reactions cannot all hold.
        KeyError
            If `stiffness` is a dict without a value for every edge.

        """
        nodes = list(self.form.nodes())
        index = {node: i for i, node in enumerate(nodes)}
        edges = list(self.form.edges())
        position = {edge: i for i, edge in enumerate(edges)}
        position.update({(v, u): i for i, (u, v) in enumerate(edges)})
        if isinstance(stiffness, dict):
            ea = [stiffness[edge] if edge in stiffness else stiffness[edge[::-1]] for edge in edges]
        else:
            ea = None if stiffness is None else [float(stiffness)] * len(edges)
        result = truss_equilibrium(
            [self.form.node_coordinates(node) for node in nodes],
            [(index[u], index[v]) for u, v in edges],
            {index[node]: list(self.form.node_load(node)) for node in nodes},
            {index[node]: self.form.node_attribute(node, "support") for node in self.form.supports()},
            stiffness=ea,
            forces={position[edge]: value for edge, value in (forces or {}).items()},
            reactions={(index[node], axis): value for (node, axis), value in (reactions or {}).items()},
            tol=tol,
        )
        for edge, value in zip(edges, result.forces):
            self.form.edge_attribute(edge, "_f", value)
        for node in nodes:
            self.form.node_attributes(node, ["_rx", "_ry", "_rz"], result.reactions.get(index[node], [0.0, 0.0, 0.0]))
        self.result = result
        self.self_stresses = [
            ({edge: force for edge, force in zip(edges, state_forces)}, {nodes[i]: Vector(*reaction) for i, reaction in state_reactions.items()})
            for state_forces, state_reactions in result.self_stresses
        ]
        self.cells, self.cell_faces, self.skipped = {}, {}, {}
        return result

    def derive_force_diagram(self, kmax: int = 2000, tolerance_flat: float = 0.001, tolerance_area: float = 0.0001) -> dict[int, Cell]:
        """Construct the force polyhedron of every node from the forces of the form diagram.

        Each force on a node, from a member, the load or the reaction, gives one face of its cell; see `cell_from_forces`.
        Nodes with fewer than four forces, or with all forces in one plane, have no polyhedron:
        they are listed in `skipped`, with the reason.

        Parameters
        ----------
        kmax
            The number of planarisation iterations per cell.
        tolerance_flat, tolerance_area
            The convergence thresholds of `cell_planarise`.

        Returns
        -------
        dict[int, Cell]
            The cell of each node.

        """
        edges = list(self.form.edges())
        key = {edge: i for i, edge in enumerate(edges)}
        load, reaction = len(edges), len(edges) + 1  # the face keys of a cell are integers
        self.cells, self.cell_faces, self.skipped = {}, {}, {}
        for node in self.form.nodes():
            forces = {}
            labels = {}
            for name, vector in self.form.node_forces(node).items():
                if TOL.is_zero(vector.length):
                    continue
                face = load if name == "load" else reaction if name == "reaction" else key[name]
                forces[face] = list(vector)
                labels[face] = name
            try:
                cell = cell_from_forces(forces, self.form.node_coordinates(node), kmax=kmax, tolerance_flat=tolerance_flat, tolerance_area=tolerance_area)
            except ValueError as error:
                self.skipped[node] = str(error)
                continue
            self.cells[node] = cell
            self.cell_faces[node] = labels
        return self.cells

    # --------------------------------------------------------------------------
    # results
    # --------------------------------------------------------------------------

    def forces(self) -> dict[tuple[int, int], float]:
        """Return the axial force of every member, tension positive.

        Returns
        -------
        dict[tuple[int, int], float]

        """
        return {edge: self.form.edge_force(edge) for edge in self.form.edges()}

    def reactions(self) -> dict[int, Vector]:
        """Return the reaction of every support.

        Returns
        -------
        dict[int, Vector]

        """
        return {node: self.form.node_reaction(node) for node in self.form.supports()}

    def check(self) -> dict[int, dict[str, float]]:
        """Measure how well each force polyhedron matches the forces on its node.

        Returns
        -------
        dict[int, dict[str, float]]
            For each node with a cell, the largest relative difference between a face area and its force (`"area_error"`),
            and the largest angle in degrees between a face normal and its force (`"angle"`).

        """
        report = {}
        for node, cell in self.cells.items():
            forces = self.form.node_forces(node)
            errors = []
            angles = []
            for face, name in self.cell_faces[node].items():
                force = forces[name]
                errors.append(abs(cell.face_area(face) - force.length) / force.length)
                angles.append(cell.face_normal(face).angle(force, degrees=True))
            report[node] = {"area_error": max(errors), "angle": max(angles)}
        return report

    def summary(self) -> str:
        """Write the member forces and support reactions as a text table.

        Returns
        -------
        str

        """
        lines = ["members (+ tension, - compression):"]
        for (u, v), force in self.forces().items():
            state = "0" if TOL.is_zero(force) else "T" if force > 0 else "C"
            lines.append("  {:>3}-{:<3} {:+10.3f}  {}".format(u, v, force, state))
        reactions = self.reactions()
        if reactions:
            lines.append("reactions:")
            for node, vector in reactions.items():
                lines.append("  node {:>3} ({:>3}): [{:+.3f}, {:+.3f}, {:+.3f}]".format(node, self.form.node_attribute(node, "support"), *vector))
        return "\n".join(lines)
