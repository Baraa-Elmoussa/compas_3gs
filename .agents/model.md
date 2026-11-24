# Model: trusses, supports, forces, equilibrium, force polyhedra

`from compas_3gs.model import GSModel` — the high-level API. It works on plain `compas_3gs`
objects (`model.form`, `model.force`, `model.cells`), so every lower-level function still applies.
Examples: `examples/model/01`–`08`.

## Form first

```python
model = GSModel.from_lines(lines)                     # [(start, end), ...]; ends joined where they meet
node = model.node_at([0.0, 0.0, 3.0])                 # KeyError if no node there (tol: TOL.absolute)
model.set_support(node=node, dofs="pinned")           # "pinned" xyz, "slider" xz, "roller" z, "free", or axes "yz"
model.add_force(node=node, vector=[0.0, 0.0, -10.0])  # adds to the force already there
result = model.solve()                                # member forces and reactions, stored in model.form
model.derive_force_diagram()                          # model.cells: one force polyhedron per node
print(model.summary())                                # text table of forces and reactions
```

| Call | Returns / effect |
|---|---|
| `GSModel(form=None, force=None, name=None)` | a `compas.data.Data`: saves to JSON with its form (and force) diagram |
| `GSModel.from_lines(lines, precision=None)` | model with a `FormNetwork` |
| `GSModel.from_force_diagram(force, kmax=500, edge_min=0.1, edge_max=1e6, tolerance=0.001)` | force first, see below |
| `model.solve(stiffness=None, forces=None, reactions=None, tol=None)` | `EquilibriumResult`; writes edge `_f`, node `_rx _ry _rz`; clears old cells |
| `model.self_stresses` | after `solve`: `[({edge: force}, {node: Vector}), ...]`, one per self-stress state |
| `model.derive_force_diagram(kmax=2000, tolerance_flat=0.001, tolerance_area=1e-4, max_area_error=0.01, max_angle=1.0)` | `{node: Cell}`; nodes without a polyhedron in `model.skipped` with the reason |
| `model.cell_faces[node]` | `{face: edge or "load" or "reaction"}`: what each cell face stands for |
| `model.check()` | `{node: {"area_error", "angle"}}`: how well each cell matches its forces |
| `model.forces()`, `model.reactions()` | `{edge: float}` (tension positive), `{node: Vector}` |

## Statically indeterminate trusses

Statics gives a family of states: one state plus any combination of the self-stress states.
`solve` picks the state with the least complementary energy `sum(N**2 L / (2 EA))` among those with
the prescribed values:

```python
model.solve()                                          # elastic state, equal EA for all members
model.solve(stiffness={edge: EA, ...})                 # elastic state for these sections (= linear FEM, rigid supports)
model.solve(forces={edge: -10.0})                      # one prescribed value per self-stress state fixes the state
model.solve(reactions={(node, "z"): 25.0})             # or prescribe reaction components
```

- `result.determinate` is True only with no self-stress state and no mechanism.
- A mechanism that cannot carry the load raises `ValueError("... mechanism ...")`; one that can
  carries it (`determinate` False, no self-stress).
- Prescribed values that statics fixes, or too many, raise `ValueError("... cannot all hold ...")`.
- Pure prestress (no loads, e.g. a tensegrity): prescribe one member force. See `08_tensegrity_prism.py`.

## Force first

```python
force = ForceVolMesh.from_vertices_and_cells(vertices, cells)
model = GSModel.from_force_diagram(force)  # form = dual network, reciprocated; forces = face areas
```

Member forces come from the face areas (`utilities.get_force_mags`), the loads from the boundary
faces. `model.form.dual is model.force` and `model.force.primal is model.form`.

## The form diagram underneath: `compas_3gs.diagrams.FormNetwork`

A `compas.datastructures.Graph`. Node attributes: `support` (axes string), `px py pz` (applied
force), `_rx _ry _rz` (reaction); edge attribute `_f` (axial force, tension positive).

| Method | |
|---|---|
| `node_at(point, tol=None)` | node at a point, or `KeyError` |
| `set_support(node, dofs="xyz")`, `supports()` | support types as above; iterate support nodes |
| `add_force(node, vector)`, `node_load(node)` | applied force, a `Vector` |
| `node_reaction(node)`, `edge_force(edge)` | after solving |
| `node_forces(node)` | `{edge: Vector, "load": Vector, "reaction": Vector}`: all forces on the node; a tension member pulls it towards its neighbour |
| `dual_halfface(edge)`, `update_angle_deviations()` | with a dual `ForceVolMesh`: the face of an edge; angle between edge and face normal into edge attribute `_a` |

## The solver on plain data: `compas_3gs.algorithms`

```python
from compas_3gs.algorithms import truss_equilibrium
result = truss_equilibrium(vertices, edges, loads, supports, stiffness=None, forces=None, reactions=None, tol=None)
```

- `vertices`: coordinates; `edges`: index pairs; `loads`: `{vertex: [fx, fy, fz]}`;
  `supports`: `{vertex: "xyz"}`; `stiffness`: EA per edge; `forces`: `{edge index: N}`;
  `reactions`: `{(vertex, "z"): R}`.
- `EquilibriumResult(forces, reactions, residual, determinate, self_stresses)`: forces per edge,
  reactions `{vertex: [rx, ry, rz]}`, self-stress states as `(forces, reactions)` scaled to a largest
  member force of 1.

`cell_from_forces(forces, origin, kmax=2000, tolerance_flat=0.001, tolerance_area=1e-4)` builds one
force polyhedron from `{key: force vector}`; it raises `ValueError` for fewer than four forces, two in
the same direction, or all in one plane.

## Limits

- Trusses only: pin-jointed members, axial forces, no moments.
- The force polyhedron of a node takes its faces and their neighbours from the force directions
  (the EGI). For nodes with several forces in one plane and of very different sizes, the polyhedron
  with those face areas has other neighbours, and the cell is not kept (`model.skipped`).
  Four forces always work (a tetrahedron).
- A force-first `GSModel` saved to JSON comes back with extra halffaces (compas `VolMesh` bug,
  [`compas2.md`](compas2.md)).
