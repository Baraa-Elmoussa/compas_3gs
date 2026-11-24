"""A compas_viewer `Viewer` for 3D graphic statics models.

This module needs `compas_viewer`, which `compas_3gs` itself does not depend on.

"""

from collections.abc import Sequence
from typing import Any

from compas_viewer import Viewer
from compas_viewer.scene import Tag

from compas.colors import Color
from compas.datastructures import Mesh
from compas.geometry import Frame
from compas.geometry import Line
from compas.geometry import Point
from compas.geometry import Scale
from compas.geometry import Translation
from compas.geometry import Vector
from compas.geometry import add_vectors
from compas.geometry import bounding_box
from compas.geometry import centroid_points
from compas.geometry import cross_vectors
from compas.geometry import distance_point_point
from compas.geometry import normalize_vector
from compas.geometry import scale_vector
from compas.geometry import subtract_vectors
from compas.tolerance import TOL
from compas_3gs.diagrams import ForceVolMesh
from compas_3gs.diagrams import FormNetwork
from compas_3gs.model import GSModel

__all__ = ["GSViewer"]

COMPRESSION = Color.red()
TENSION = Color.blue()
ZERO = Color.grey()
LOAD = Color.black()
REACTION = Color.from_hex("#2e9e44")
SUPPORT = Color.from_hex("#6d6d6d")
AXIS = {"x": [1.0, 0.0, 0.0], "y": [0.0, 1.0, 0.0], "z": [0.0, 0.0, 1.0]}


class GSViewer(Viewer):
    """A viewer that draws a `GSModel`: the form diagram with its supports, loads and reactions, and the force diagram.

    Everything goes into groups of the scene tree, so that each part can be shown, hidden or removed on its own:

    - `form`: `members` (red in compression, blue in tension, the width by the force), `labels` (the member forces),
      `supports` (a glyph per support), `loads` (black arrows) and `reactions` (green arrows);
    - `force`: `cells`, the force polyhedron of every node, or the polyhedral force diagram.

    A support is drawn as a tetrahedron under the node, on two rollers if it restrains one direction,
    and on a line along the free direction if it restrains two.

    Examples
    --------
    ```python
    viewer = GSViewer()
    viewer.add_model(model)
    viewer.show()
    ```

    """

    def add_model(self, model: GSModel, labels: bool = True, offset: Sequence[float] | None = None) -> dict[str, Any]:
        """Add a model: its form diagram and its force diagram or force polyhedra.

        Parameters
        ----------
        model
            The model.
        labels
            Show the member forces as labels.
        offset
            Where the force diagram goes relative to the form diagram. Default is beside it, along x.

        Returns
        -------
        dict[str, Any]
            The groups, by name: `"form"`, `"members"`, `"labels"`, `"supports"`, `"loads"`, `"reactions"`
            and, if the model has a force diagram or cells, `"force"` and `"cells"`.

        """
        groups = self.add_form(model.form, labels=labels)
        if offset is None:
            box = bounding_box([model.form.node_coordinates(node) for node in model.form.nodes()])
            offset = [1.5 * (box[6][0] - box[0][0]) + 2.0, 0.0, 0.0]
        if model.force is not None:
            groups.update(self.add_force_diagram(model.force, model.form, offset))
        elif model.cells:
            groups.update(self.add_cells(model, offset))
        return groups

    def add_form(self, form: FormNetwork, labels: bool = True) -> dict[str, Any]:
        """Add a form diagram with its members, supports, loads and reactions.

        Parameters
        ----------
        form
            The form diagram.
        labels
            Show the member forces as labels.

        Returns
        -------
        dict[str, Any]
            The groups, by name: `"form"`, `"members"`, `"labels"`, `"supports"`, `"loads"` and `"reactions"`.

        """
        box = bounding_box([form.node_coordinates(node) for node in form.nodes()])
        size = 0.05 * distance_point_point(box[0], box[6]) or 1.0

        group = self.scene.add_group(name="form")
        members = self.scene.add_group(name="members", parent=group)
        texts = self.scene.add_group(name="labels", parent=group)
        supports = self.scene.add_group(name="supports", parent=group)
        loads = self.scene.add_group(name="loads", parent=group)
        reactions = self.scene.add_group(name="reactions", parent=group)

        forces = {edge: form.edge_force(edge) for edge in form.edges()}
        peak = max([abs(force) for force in forces.values()], default=0.0) or 1.0
        for (u, v), force in forces.items():
            a, b = form.node_coordinates(u), form.node_coordinates(v)
            color = ZERO if TOL.is_zero(force) else TENSION if force > 0 else COMPRESSION
            self.scene.add(Line(a, b), parent=members, linecolor=color, linewidth=1.5 + 6 * abs(force) / peak, name="{}-{}".format(u, v))
            if labels:
                self.scene.add(Tag("{:+.2f}".format(force), centroid_points([a, b]), color=color, height=14), parent=texts)

        lengths = [form.node_load(node).length for node in form.nodes()] + [form.node_reaction(node).length for node in form.supports()]
        arrow = 4 * size / (max(lengths, default=0.0) or 1.0)  # the longest arrow is four glyph sizes long
        for node in form.nodes():
            point = form.node_coordinates(node)
            load = form.node_load(node)
            if not TOL.is_zero(load.length):
                vector = Vector(*scale_vector(load, arrow))
                self.scene.add(vector, parent=loads, anchor=Point(*subtract_vectors(point, vector)), linecolor=LOAD, linewidth=2)
            dofs = form.node_attribute(node, "support")
            if dofs:
                self.add_support(point, dofs, size, parent=supports)
                reaction = form.node_reaction(node)
                if not TOL.is_zero(reaction.length):
                    vector = Vector(*scale_vector(reaction, arrow))
                    self.scene.add(vector, parent=reactions, anchor=Point(*subtract_vectors(point, vector)), linecolor=REACTION, linewidth=2)
        return {"form": group, "members": members, "labels": texts, "supports": supports, "loads": loads, "reactions": reactions}

    def add_support(self, point: Sequence[float], dofs: str, size: float, parent: Any = None) -> None:
        """Add the glyph of a support: a tetrahedron under the point, along its main restrained direction.

        Parameters
        ----------
        point
            The location of the support.
        dofs
            The restrained directions, as axes, for example `"xyz"`.
        size
            The height of the tetrahedron.
        parent
            The group of the glyph.

        """
        axis = "z" if "z" in dofs else dofs[0]
        d = AXIS[axis]
        e1 = normalize_vector(cross_vectors(d, AXIS["y"] if axis != "y" else AXIS["x"]))
        e2 = cross_vectors(d, e1)
        base = subtract_vectors(point, scale_vector(d, size))
        corners = [add_vectors(base, add_vectors(scale_vector(e1, 0.7 * size * c), scale_vector(e2, 0.7 * size * s))) for c, s in ((1, 0), (-0.5, 0.866), (-0.5, -0.866))]
        tetrahedron = Mesh.from_vertices_and_faces([point] + corners, [[0, 1, 2], [0, 2, 3], [0, 3, 1], [1, 3, 2]])
        self.scene.add(tetrahedron, parent=parent, facecolor=SUPPORT, linecolor=Color.black(), show_lines=True)
        below = subtract_vectors(base, scale_vector(d, 0.25 * size))
        if len(dofs) == 1:  # a roller: two rollers under the base
            half = scale_vector(e2, 0.5 * size)
            for offset in (-0.4, 0.4):
                centre = add_vectors(below, scale_vector(e1, offset * size))
                self.scene.add(Line(subtract_vectors(centre, half), add_vectors(centre, half)), parent=parent, linecolor=Color.black(), linewidth=3)
        elif len(dofs) == 2:  # a slider: a line along the free direction
            free = scale_vector(AXIS[next(a for a in "xyz" if a not in dofs)], size)
            self.scene.add(Line(subtract_vectors(below, free), add_vectors(below, free)), parent=parent, linecolor=Color.black(), linewidth=3)

    def add_cells(self, model: GSModel, offset: Sequence[float]) -> dict[str, Any]:
        """Add the force polyhedron of every node, at the node moved by `offset`.

        The cells are scaled together, so that they do not overlap. A face has the colour of its member,
        black for a load, green for a reaction and white if it carries no force.

        Parameters
        ----------
        model
            The model, after `derive_force_diagram`.
        offset
            The translation from the form diagram to the cells.

        Returns
        -------
        dict[str, Any]
            The groups `"force"` and `"cells"`.

        """
        form = model.form
        group = self.scene.add_group(name="force")
        cells = self.scene.add_group(name="cells", parent=group)
        diameters = []
        for cell in model.cells.values():
            points = [cell.vertex_coordinates(vertex) for vertex in cell.vertices()]
            diameters.append(max(distance_point_point(a, b) for a in points for b in points))
        shortest = min(form.edge_length(edge) for edge in form.edges())
        k = 0.6 * shortest / (max(diameters, default=0.0) or 1.0)
        for node, cell in model.cells.items():
            centre = centroid_points([cell.vertex_coordinates(vertex) for vertex in cell.vertices()])
            target = add_vectors(form.node_coordinates(node), offset)
            T = Translation.from_vector(subtract_vectors(target, centre)) * Scale.from_factors([k, k, k], frame=Frame(centre, [1, 0, 0], [0, 1, 0]))
            facecolor = {}
            for face in cell.faces():
                name = model.cell_faces[node].get(face)
                if name is None:
                    facecolor[face] = Color.white()
                elif name == "load":
                    facecolor[face] = LOAD
                elif name == "reaction":
                    facecolor[face] = REACTION
                else:
                    force = form.edge_force(name)
                    facecolor[face] = ZERO if TOL.is_zero(force) else TENSION if force > 0 else COMPRESSION
            self.scene.add(cell.transformed(T), parent=cells, facecolor=facecolor, opacity=0.75, show_lines=True, name="cell {}".format(node))
        return {"force": group, "cells": cells}

    def add_force_diagram(self, force: ForceVolMesh, form: FormNetwork, offset: Sequence[float]) -> dict[str, Any]:
        """Add a polyhedral force diagram, one mesh per cell, its faces coloured like the members of the form diagram.

        Boundary faces, the loads, are black.

        Parameters
        ----------
        force
            The force diagram.
        form
            Its form diagram, with the member forces.
        offset
            The translation from the centre of the form diagram to the centre of the force diagram.

        Returns
        -------
        dict[str, Any]
            The groups `"force"` and `"cells"`.

        """
        group = self.scene.add_group(name="force")
        cells = self.scene.add_group(name="cells", parent=group)
        centre = centroid_points([form.node_coordinates(node) for node in form.nodes()])
        T = Translation.from_vector(subtract_vectors(add_vectors(centre, offset), force.centroid()))
        colors = {}
        for u, v in form.edges():
            halfface = force.cell_pair_halffaces(u, v)[0]
            value = form.edge_force((u, v))
            colors[frozenset(force.halfface_vertices(halfface))] = ZERO if TOL.is_zero(value) else TENSION if value > 0 else COMPRESSION
        for cell in force.cells():
            mesh = force.cell_to_mesh(cell)
            vertices = force.cell_vertices(cell)
            facecolor = {face: colors.get(frozenset(vertices[i] for i in mesh.face_vertices(face)), LOAD) for face in mesh.faces()}
            self.scene.add(mesh.transformed(T), parent=cells, facecolor=facecolor, opacity=0.6, show_lines=True, name="cell {}".format(cell))
        return {"force": group, "cells": cells}
