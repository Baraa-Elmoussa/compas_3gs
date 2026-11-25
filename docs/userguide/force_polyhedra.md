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

By Minkowski's theorem, face normals that span space and face areas whose vectors (area times normal) sum to
zero belong to exactly one convex polyhedron. Forces in equilibrium are such faces, so every node has one
force polyhedron, and which faces of it meet follows from the force sizes as much as from their directions.

`cell_from_forces` finds it as Minkowski proved it exists. The distances of the face planes from the centre that
give the largest volume for a given weighted sum of the distances also give face areas in proportion to the
forces; a quasi-Newton optimisation finds them, the cell is scaled, and Newton steps bring every face to the
area of its force. For one node:

```python
from compas_3gs.algorithms import cell_from_forces, cell_from_normals_and_areas

cell = cell_from_forces({0: f0, 1: f1, 2: f2, 3: f3}, origin)  # the keys become the face keys
cell = cell_from_normals_and_areas(normals, areas)              # the same, from normals and areas
```

The [polyhedral cell constructor](polyhedral_cell/00_20_cell_constructor.md) of the user manual builds a cell
in another way: the force directions give the extended Gaussian image (EGI), the EGI a cell with those normals,
and planarisation sizes its faces ([Planarisation](multi_cell_polyhedron/01_40_planarisation.md)). Its faces
and their neighbours come from the directions alone, so it suits cells whose faces are designed, and the
functions it uses stay in `compas_3gs.algorithms`.

## Which nodes have a polyhedron

- A node needs four forces or more, not all in one plane, and in equilibrium.
- No two of its forces may point the same way: a convex polyhedron cannot have two faces with one normal. Two
  collinear members at a node, one pulling and one pushing in the same direction, are such a pair
  ([Transmission tower](../examples/model_06_transmission_tower.md): the legs run straight through the levels).
- Four forces give a tetrahedron ([Tensegrity prism](../examples/model_08_tensegrity_prism.md)); a lattice node
  with nine forces of very different sizes gives a thin slab with small faces, and both are exact.

## Force first

A model built from a polyhedral force diagram already has its force diagram: `model.force`, one cell per
node, the cells sharing the faces of the members. Its cells are the subject of the
[multi-cell polyhedron](multi_cell_polyhedron/index.md) pages.

Next: [Viewing](viewing.md).
