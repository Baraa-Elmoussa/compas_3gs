# Examples

These scripts show how to use the API of `compas_3gs`, with `compas_viewer` for the drawings.
Run them from the root of the repository, for example `python examples/model/01_tripod.py`.

## The API in brief

A `GSModel` holds a form diagram with its supports and forces, solves it, and builds its force diagram;
a `GSViewer` draws all of it, grouped. The diagrams stay plain `compas_3gs` objects (`model.form`, `model.force`,
`model.cells`), so every algorithm and operation of the package works on them, as the other examples show.

```python
from compas_3gs.model import GSModel
from compas_3gs.viewer import GSViewer

# form first: member lines, supports and forces
model = GSModel.from_lines(lines)
model.set_support(node=model.node_at(point), dofs="pinned")  # or "slider", "roller", or axes such as "xz"
model.add_force(node=model.node_at(apex), vector=[0.0, 0.0, -10.0])
model.solve()  # member forces and support reactions from equilibrium
model.derive_force_diagram()  # one force polyhedron per node
print(model.summary())

# force first: a polyhedral force diagram gives its form diagram
model = GSModel.from_force_diagram(volmesh)

viewer = GSViewer()
viewer.add_model(model)
viewer.show()
```

## Models

A form diagram from lines, supports and forces, solved with `GSModel` and drawn with `GSViewer`; and a model from a force diagram.

- [A tripod with `GSModel`: member lines, supports and a load, solve, force polyhedra, view](model_01_tripod.md)
- [A tetrahelix tower on a pinned, a slider and a roller support, with `GSModel`](model_02_space_truss.md)
- [Force first with `GSModel`: a polyhedral force diagram gives the frame and its forces](model_03_force_first.md)
- [What `GSModel.solve` does when statics alone has many solutions, or none](model_04_determinacy.md)
- [A 3D Warren truss solved with 3GS (`GSModel`, statics) and with a linear finite element model](model_05_warren_truss_fem.md)

## Physical structures

The same kind of structures built step by step with the algorithms: equilibrium, EGI, cells, duality and reciprocation.

- [A tripod: three legs carrying a load at the apex, and its force polyhedron](physical_structures_10_tripod.md)
- [A space truss tower (a tetrahelix) and one force polyhedron per node](physical_structures_20_space_truss.md)
- [Force first: design a polyhedral force diagram, get the space frame it stands for](physical_structures_30_force_first_frames.md)
- [Thrust line of a semicircular masonry arch, and what 3GS can and cannot do with it](physical_structures_40_arch_thrust_line.md)

## Polyhedral cell

Constructing and manipulating a single polyhedral cell, the force polyhedron of one node.

- [Pull a face of a cell along its normal](cell_manipulation_00_10_cell_pull_face.md)
- [Arearise one face of a cell: move it along its normal until it reaches a target area](cell_manipulation_00_20_cell_arearise_face.md)
- [Construct a cell from force vectors in equilibrium](cell_manipulation_00_30_cell_constructor.md)
- [Planarise the faces of a cell](cell_manipulation_00_40_cell_planarisation.md)
- [Arearise several faces of a cell at once, optionally keeping all face normals fixed](cell_manipulation_00_41_cell_arearisation.md)

## Multi-cell polyhedron

Polyhedral force diagrams of several cells, their dual form diagrams, reciprocation, planarisation and visualisation.

- [Construct a force volmesh from a set of cells and show its cells, vertices and vertex keys](multi_cell_polyhedron_01_10_volmesh_constructor.md)
- [Construct the dual volmesh of a force volmesh](multi_cell_polyhedron_01_20_volmesh_duality_volmesh.md)
- [Make the dual network of a force volmesh and show how their elements correspond](multi_cell_polyhedron_01_21_volmesh_duality_network.md)
- [Reciprocate a polyhedral force diagram and its form diagram](multi_cell_polyhedron_01_30_volmesh_reciprocation.md)
- [Planarise the faces of a volmesh, with some of its vertices fixed](multi_cell_polyhedron_01_40_volmesh_planarisation.md)
- [Arearise chosen halffaces of a volmesh to a target area](multi_cell_polyhedron_01_50_volmesh_arearisation.md)
- [Modify a polyhedral force diagram: set vertex fixities, move a vertex and lift a boundary vertex](multi_cell_polyhedron_01_60_volmesh_modification.md)
- [Colour the corresponding elements of a polyhedral form and force diagram](multi_cell_polyhedron_01_70_volmesh_color_id.md)
- [Five ways to visualise a polyhedral force diagram and its reciprocal form diagram](multi_cell_polyhedron_01_71_volmesh_visualisation.md)
- [Unified diagram of a multi-cell force diagram and its form diagram](multi_cell_polyhedron_01_80_volmesh_unified_diagram.md)
