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

    # the data methods, `delete_cell`, `is_valid`, `edge_halffaces`, `halfface_manifold_neighbors`,
    # `cell_vertex_neighbors` and `cell_vertex_faces` replace VolMesh methods that fail in compas 2.15 and 2.16 (TODO.md)

    @property
    def __data__(self) -> dict[str, Any]:
        cells = {cell: sorted(self.cell_faces(cell)) for cell in self.cells()}
        return {
            "attributes": self.attributes,
            "default_vertex_attributes": self.default_vertex_attributes,
            "default_edge_attributes": self.default_edge_attributes,
            "default_face_attributes": self.default_face_attributes,
            "default_cell_attributes": self.default_cell_attributes,
            "vertex": {str(vertex): attr for vertex, attr in self._vertex.items()},
            "cell": {str(cell): [self.halfface_vertices(halfface) for halfface in halffaces] for cell, halffaces in cells.items()},
            "cell_halffaces": {str(cell): halffaces for cell, halffaces in cells.items()},
            "edge_data": self._edge_data,
            "face_data": self._face_data,
            "cell_data": {str(cell): attr for cell, attr in self._cell_data.items()},
            "max_vertex": self._max_vertex,
            "max_face": self._max_face,
            "max_cell": self._max_cell,
        }

    @classmethod
    def __from_data__(cls, data: dict[str, Any]) -> "VolMesh3gs":
        volmesh = cls(
            default_vertex_attributes=data.get("default_vertex_attributes"),
            default_edge_attributes=data.get("default_edge_attributes"),
            default_face_attributes=data.get("default_face_attributes"),
            default_cell_attributes=data.get("default_cell_attributes"),
        )
        volmesh.attributes.update(data["attributes"] or {})

        for vertex, attr in (data["vertex"] or {}).items():
            volmesh.add_vertex(key=vertex, attr_dict=attr)

        cell_data = data.get("cell_data") or {}
        cell_halffaces = data.get("cell_halffaces") or {}
        for key, faces in (data["cell"] or {}).items():
            # add_cell, with the halfface keys of the data
            cell = volmesh.add_cell([], ckey=key, attr_dict=cell_data.get(key))
            added = set()
            for halfface, vertices in zip(cell_halffaces.get(key) or [None] * len(faces), faces):
                if frozenset(vertices) in added:
                    continue  # the data of a compas VolMesh lists every face once per halfedge
                added.add(frozenset(vertices))
                halfface = volmesh.add_halfface(vertices, fkey=halfface)
                vertices = volmesh.halfface_vertices(halfface)
                for i in range(-2, len(vertices) - 2):
                    u, v, w = vertices[i], vertices[i + 1], vertices[i + 2]
                    volmesh._cell[cell].setdefault(u, {})[v] = halfface
                    volmesh._plane[u][v][w] = cell

        for edge, attr in (data.get("edge_data") or {}).items():
            volmesh._edge_data[edge] = attr or {}
        for face, attr in (data.get("face_data") or {}).items():
            volmesh._face_data[face] = attr or {}

        volmesh._max_vertex = data.get("max_vertex", volmesh._max_vertex)
        volmesh._max_face = data.get("max_face", volmesh._max_face)
        volmesh._max_cell = data.get("max_cell", volmesh._max_cell)
        return volmesh

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
    #   topology
    # --------------------------------------------------------------------------

    def delete_cell(self, cell: int) -> None:
        """Delete a cell, with the attributes of the faces and edges that no other cell uses.

        Parameters
        ----------
        cell
            The identifier of the cell.

        Raises
        ------
        KeyError
            If the cell does not exist.

        """
        faces = [str(tuple(sorted(self.halfface_vertices(halfface)))) for halfface in self.cell_faces(cell) if self.halfface_opposite_cell(halfface) is None]
        edges = [str(tuple(sorted(edge))) for edge in self.cell_edges(cell) if set(self.edge_cells(edge)) == {cell}]
        super().delete_cell(cell)
        for face in faces:
            self._face_data.pop(face, None)
        for edge in edges:
            self._edge_data.pop(edge, None)

    def is_valid(self) -> bool:
        """Verify that every cell is closed and that neighbouring cells meet in matching faces.

        Returns
        -------
        bool
            True if every halfface belongs to a cell and has at least three vertices, every halfedge of a cell has its
            opposite in the same cell, and every halfface shared with a neighbouring cell has the same vertices as the
            halfface of that cell, in the opposite order.

        """
        for halfface, vertices in self._halfface.items():
            if len(set(vertices)) < 3:
                return False
            u, v = vertices[:2]
            cell = self.halfface_cell(halfface)
            if cell is None or self._cell[cell].get(u, {}).get(v) != halfface:
                return False
            nbr = self.halfface_opposite_cell(halfface)
            if nbr is None:
                continue
            opposite = self._cell[nbr].get(v, {}).get(u)
            if opposite is None:
                return False
            cycle = self.halfface_vertices(opposite)[::-1]
            if len(cycle) != len(vertices) or u not in cycle:
                return False
            i = cycle.index(u)
            if cycle[i:] + cycle[:i] != vertices:
                return False

        for cell in self.cells():
            for u in self._cell[cell]:
                for v in self._cell[cell][u]:
                    if u not in self._cell[cell].get(v, {}):
                        return False
        return True

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
    #   edges and halffaces
    # --------------------------------------------------------------------------

    def edge_halffaces(self, edge: tuple[int, int]) -> list[int]:
        """Return the halffaces running along an edge, one per cell, in order around the edge.

        Parameters
        ----------
        edge
            The identifier of the edge, `(u, v)`.

        Returns
        -------
        list[int]
            For every cell around the edge, its halfface with the halfedge `(u, v)`. The next cell is the one on the other
            side of that halfface; around an edge on the boundary, the first cell is on the boundary.

        """
        u, v = edge
        cells = [cell for cell in self._plane[u][v].values() if cell is not None]
        halffaces = []
        while cells:
            # step back to the first cell of the fan, through the halffaces with the halfedge (v, u)
            cell = cells[0]
            for _ in range(len(cells) - 1):
                previous = self.halfface_opposite_cell(self._cell[cell][v][u])
                if previous is None or previous == cells[0]:
                    break
                cell = previous
            while cell in cells:
                halfface = self._cell[cell][u][v]
                halffaces.append(halfface)
                cells.remove(cell)
                cell = self.halfface_opposite_cell(halfface)
        return halffaces

    def halfface_manifold_neighbors(self, halfface: int) -> list[int]:
        """Return the halffaces that continue a halfface across its edges, in the neighbouring cells.

        Parameters
        ----------
        halfface
            The identifier of the halfface.

        Returns
        -------
        list[int]
            For every edge of the halfface where the cell continues into a neighbouring cell, the halfface of that
            neighbour along the edge. Halffaces of the same cell are not included.

        """
        cell = self.halfface_cell(halfface)
        nbrs = []
        for u, v in self.halfface_halfedges(halfface):
            nbr = self.halfface_opposite_cell(self._cell[cell][v][u])
            if nbr is not None:
                nbrs.append(self._cell[nbr][v][u])
        return nbrs

    # --------------------------------------------------------------------------
    # cell
    # --------------------------------------------------------------------------

    def cell_vertex_neighbors(self, cell: int, vertex: int) -> list[int]:
        """Return the neighbours of a vertex in a cell, in order around the vertex.

        Parameters
        ----------
        cell
            The identifier of the cell.
        vertex
            The identifier of the vertex.

        Returns
        -------
        list[int]
            The vertices connected to the vertex by an edge of the cell.

        Raises
        ------
        KeyError
            If the vertex is not a vertex of the cell.

        """
        if vertex not in self._cell[cell]:
            raise KeyError(vertex)
        nbr = next(iter(self._cell[cell][vertex]))
        nbrs = [nbr]
        for _ in range(len(self._cell[cell][vertex]) - 1):
            nbr = self.halfface_vertex_ancestor(self._cell[cell][vertex][nbr], vertex)
            nbrs.append(nbr)
        return nbrs

    def cell_vertex_faces(self, cell: int, vertex: int) -> list[int]:
        """Return the halffaces of a cell around one of its vertices, in order around the vertex.

        Parameters
        ----------
        cell
            The identifier of the cell.
        vertex
            The identifier of the vertex.

        Returns
        -------
        list[int]
            The halffaces of the cell with the vertex; halfface `i` has the halfedge from the vertex to neighbour `i` of
            `cell_vertex_neighbors(cell, vertex)`.

        Raises
        ------
        KeyError
            If the vertex is not a vertex of the cell.

        """
        return [self._cell[cell][vertex][nbr] for nbr in self.cell_vertex_neighbors(cell, vertex)]

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
