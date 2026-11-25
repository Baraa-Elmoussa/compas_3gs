# Construct the dual volmesh of a force volmesh

Every cell of the force volmesh becomes a vertex of the dual volmesh (same key, same colour), and
every interior vertex a cell. Each dual edge (u, v) is drawn as an arrow from u to v, and the halfface
of cell u that points to cell v as a face with its normal, in the same colour as the edge.

![Construct the dual volmesh of a force volmesh](../assets/images/examples/multi_cell_polyhedron_01_20_volmesh_duality_volmesh.png)

```python title="examples/multi_cell_polyhedron/01_20_volmesh_duality_volmesh.py"
--8<-- "examples/multi_cell_polyhedron/01_20_volmesh_duality_volmesh.py"
```
