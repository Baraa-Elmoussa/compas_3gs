# Tasks: what to call, and the example that does it

Look here before writing code: these all exist. Examples run from the repository root.

| I want to | Call | Example |
|---|---|---|
| model a truss from member lines | `GSModel.from_lines(lines)` | `examples/model/01_tripod.py` |
| add supports (pinned, slider, roller, any axes) | `model.set_support(node=model.node_at(p), dofs="pinned")` | `examples/model/02_space_truss.py` |
| add loads | `model.add_force(node=..., vector=[...])` | `examples/model/01_tripod.py` |
| get member forces and reactions | `model.solve()`, `model.forces()`, `model.reactions()`, `model.summary()` | `examples/model/01_tripod.py` |
| know if statics is enough | `result.determinate`, `len(model.self_stresses)`; mechanisms raise `ValueError` | `examples/model/04_determinacy.py` |
| solve an indeterminate truss elastically | `model.solve(stiffness={edge: EA})` | `examples/model/06_transmission_tower.py` |
| pick another equilibrium state | `model.solve(forces={edge: N})` or `reactions={(node, "z"): R}` | `examples/model/05_warren_truss_fem.py` |
| prestress a tensegrity or cable structure | `model.solve(forces={strut: -N})` with no loads | `examples/model/08_tensegrity_prism.py` |
| compare with linear finite elements | the elastic `solve` is the FEM state; an inline FEM check | `examples/model/05_warren_truss_fem.py` |
| get the force polyhedron of every node | `model.derive_force_diagram()`, `model.cells`, `model.check()`, `model.skipped` | `examples/model/02_space_truss.py` |
| get one force polyhedron from force vectors | `cell_from_forces(forces, origin)` | `examples/physical_structures/10_tripod.py` |
| go from a force diagram to its structure | `GSModel.from_force_diagram(volmesh)` | `examples/model/03_force_first.py` |
| make a form diagram reciprocal to a force diagram | `volmesh_dual_network(volmesh, cls=FormNetwork)`, `volmesh_reciprocate(volmesh, form)` | `examples/multi_cell_polyhedron/01_30_volmesh_reciprocation.py` |
| switch compression and tension | `point_reflection(volmesh)` | `examples/physical_structures/30_force_first_frames.py` |
| build a cell from an EGI | `egi_from_vectors(vectors, origin)`, `cell_from_egi(egi)`, `cell_planarise(...)` | `examples/cell_manipulation/00_30_cell_constructor.py` |
| make faces planar, or reach target areas | `cell_planarise`, `volmesh_planarise`, `cell_arearise_face` | `examples/cell_manipulation/00_40_cell_planarisation.py`, `examples/multi_cell_polyhedron/01_40_volmesh_planarisation.py` |
| pull a face of a cell | `cell_relocate_face(cell, fkey, xyz, normal)` | `examples/cell_manipulation/00_10_cell_pull_face.py` |
| edit a force diagram (move, lift, pinch, merge, subdivide) | `compas_3gs.operations` | `examples/multi_cell_polyhedron/01_60_volmesh_modification.py` |
| draw a unified diagram | `volmesh_ud(volmesh, network, scale)` | `examples/multi_cell_polyhedron/01_80_volmesh_unified_diagram.py` |
| colour form and force diagrams by force | `get_force_colors_uv`, `get_force_colors_hf` | `examples/multi_cell_polyhedron/01_70_volmesh_color_id.py` |
| draw a model, select members to see their forces | `GSViewer().add_model(model)` | `examples/model/06_transmission_tower.py` |
| update a drawing live (slider, animation) | `viewer.update_form(form, groups)`, `@viewer.on(...)` | `examples/model/05_warren_truss_fem.py` |
| load the sample force diagrams and cells | `data/*.json`, see [`diagrams.md`](diagrams.md) | `examples/multi_cell_polyhedron/01_10_volmesh_constructor.py` |

Thrust lines of arches and other problems outside pin-jointed trusses: see
`examples/physical_structures/40_arch_thrust_line.py` for what 3GS can and cannot do there.
