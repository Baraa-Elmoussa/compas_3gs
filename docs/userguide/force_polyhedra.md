# Force polyhedra

The force diagram of a solved form diagram has one polyhedral cell per node. Each force on the node,
from a member, the load or the reaction, gives one face of the cell: normal to the force, with an area
equal to the force's magnitude. The faces close into a polyhedron because the forces are in equilibrium.
{ .lead }

```python
model.solve()
cells = model.derive_force_diagram()  # {node: Cell}
```

| Attribute | |
|---|---|
| `model.cells` | the force polyhedron of each node, a `Cell` |
| `model.cell_faces[node]` | what each face of the cell stands for: a member (its edge), `"load"` or `"reaction"`; other faces carry no force |
| `model.check()` | for each cell, the largest relative difference between a face area and its force (`"area_error"`), and the largest angle between a face normal and its force (`"angle"`, degrees) |
| `model.skipped` | the nodes without a polyhedron, with the reason |

## How a cell is built

The construction is that of the [polyhedral cell constructor](polyhedral_cell/00_20_cell_constructor.md):
the force directions give the extended Gaussian image (EGI) of the cell, the EGI gives a cell with the
right face normals, and planarisation brings each face to the area of its force
([Planarisation](multi_cell_polyhedron/01_40_planarisation.md)). For one node:

```python
from compas_3gs.algorithms import cell_from_forces

cell = cell_from_forces({0: f0, 1: f1, 2: f2, 3: f3}, origin)  # the keys become the face keys
```

## Which nodes have a polyhedron

- A node needs four forces or more, not all in one plane, and no two in the same direction.
- A cell is kept only if it reaches its forces: by default within 1 % of the area and 1° of the direction
  (`derive_force_diagram(max_area_error=0.01, max_angle=1.0)`).
- The EGI fixes the faces of the cell and which faces touch from the force directions alone. For a node with
  several forces in one plane and of very different sizes, the polyhedron with those face areas has faces
  that touch differently, and planarisation cannot reach it: the node is skipped.
- Four forces always give a tetrahedron ([Tensegrity prism](../examples/model_08_tensegrity_prism.md));
  many lattice nodes are skipped ([Transmission tower](../examples/model_06_transmission_tower.md)).

## Force first

A model built from a polyhedral force diagram already has its force diagram: `model.force`, one cell per
node, the cells sharing the faces of the members. Its cells are the subject of the
[multi-cell polyhedron](multi_cell_polyhedron/index.md) pages.

Next: [Viewing](viewing.md).
