# The model

A `GSModel` holds a structure as a form diagram, its supports and loads, its equilibrium and its force
diagram. It sits on top of the classes and algorithms of the other pages of this guide; its diagrams
stay plain `compas_3gs` objects, so every algorithm and operation of the package works on them.
{ .lead }

```python
from compas_3gs.model import GSModel
```

## Form first: from member lines

```python
model = GSModel.from_lines(lines)  # [(start, end), ...], joined where the ends meet
```

`model.form` is a `FormNetwork`: one node per joint and one edge per member, pin-jointed. The members
carry axial forces only.

`model.node_at(point)` finds the node at a point (within `TOL.absolute`, or a `tol` of your own); it
raises `KeyError` when there is none.

## Supports

```python
model.set_support(node=node, dofs="pinned")
```

| `dofs` | Restrained directions | Reaction components |
|---|---|---|
| `"pinned"` | x, y, z | 3 |
| `"slider"` | x, z | 2 |
| `"roller"` | z | 1 |
| `"free"` | none | 0 |
| axes, for example `"yz"` | the axes named | as many |

The members carry no moments, so a fixed support is a pinned one. A support is stored in the node
attribute `support` as its axes (`"xyz"`, `"xz"`, ...); `model.form.supports()` iterates over the
support nodes.

## Loads

```python
model.add_force(node=node, vector=[0.0, 0.0, -10.0])
```

A force adds to the force already on the node. Loads are the node attributes `px`, `py` and `pz`;
`model.form.node_load(node)` returns them as a vector. No units are imposed; the examples use kN and m.

## Force first: from a polyhedral force diagram

```python
from compas_3gs.diagrams import ForceVolMesh

force = ForceVolMesh.from_vertices_and_cells(vertices, cells)
model = GSModel.from_force_diagram(force)
```

The form diagram is the dual network of the force diagram: one node per cell, one member per face
shared by two cells. It is reciprocated, so that every member is perpendicular to its face
([Reciprocation](multi_cell_polyhedron/01_30_reciprocation.md)). The member forces are the face areas,
and the boundary faces of each cell are the load on its node. `model.form.dual` is the force diagram
and `model.force.primal` the form diagram.

## Saving

A `GSModel` is a compas `Data` object: `model.to_json(path)` and `GSModel.from_json(path)` save and
load it with its form diagram (supports, loads and forces included) and its force diagram.

Next: [Solving](solving.md).
