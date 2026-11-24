# Algorithms

`from compas_3gs.algorithms import ...` — plain functions; most change their input in place.

## Cells (one node)

| Function | Does |
|---|---|
| `egi_from_vectors(vectordict, origin, tol=None)` | `EGI` of a set of vectors `{int key: vector}`; raises `Exception` for fewer than four, coincident or non-spanning vectors (coplanar input is not caught, see `TODO.md`) |
| `cell_from_egi(egi)` | unit `Cell` of an EGI; face keys = EGI vertex keys = the vector keys |
| `cell_from_forces(forces, origin, kmax=2000, tolerance_flat=0.001, tolerance_area=1e-4)` | force polyhedron: EGI → cell → planarise to areas = magnitudes; `ValueError` for unusable force sets |
| `cell_planarise(cell, kmax=100, target_centers=None, target_normals=None, target_areas=None, fix_vkeys=None, avg_fkeys=None, collapse_edge_length=0.1, tolerance_flat=0.001, tolerance_area=0.001, callback=None, callback_args=None, print_result_info=False)` | planar faces, optionally with target normals, centres and areas (by face key) |
| `cell_arearise_face(cell, fkey, target_area, tol=0.001, callback=None, callback_args=None)` | move one face along its normal to a target area (golden section search; see the limits in `TODO.md`) |
| `golden_section_search(function, a, b, tol=0.001)` | minimum of a unimodal function on `[a, b]` |

## Multi-cell force diagrams (volmeshes)

| Function | Does |
|---|---|
| `volmesh_dual_volmesh(volmesh, cls=None)` | dual volmesh: a vertex per cell, a cell per interior vertex |
| `volmesh_dual_network(volmesh, cls=None)` | dual network (the form diagram): a node per cell at its centroid, an edge per pair of cells that share a face; pass `cls=FormNetwork` |
| `volmesh_reciprocate(volmesh, formdiagram, kmax=100, weight=1.0, fix_vkeys=None, edge_min=0.1, edge_max=10000, tolerance=1.0, tolerance_boundary=1.0, callback=None, callback_args=None, print_result_info=False)` | make every form edge perpendicular to its face: `weight=1` moves the form diagram only, `0` the force diagram only; tolerances are angles in degrees; boundary face orientations stay fixed |
| `volmesh_planarise(volmesh, kmax=100, target_centers=None, target_normals=None, target_areas=None, fix_vkeys=None, fix_boundary_normals=False, fix_all_normals=False, min_area=None, max_area=None, tolerance_flat=0.001, tolerance_area=0.001, tolerance_perp=0.001, callback=None, callback_args=None, print_result_info=False)` | planar halffaces of a `VolMesh3gs`/`ForceVolMesh`, with targets keyed by either halfface of a face; `min_area`, `max_area` are not used |
| `volmesh_ud(volmesh, network, scale=0.5)` | unified diagram: `(cells, prisms)`, each `{key: {"vertices", "faces"}}`, cells scaled towards their nodes, one prism per edge |
| `cellnetwork_ud(cellnetwork)` | raises `NotImplementedError` |

## Equilibrium

| Function | Does |
|---|---|
| `truss_equilibrium(vertices, edges, loads, supports, stiffness=None, forces=None, reactions=None, tol=None)` | member forces and reactions; see [`model.md`](model.md) |
| `EquilibriumResult` | `(forces, reactions, residual, determinate, self_stresses)` |

Callbacks are called after every iteration as `callback(cell, k, callback_args)` (`cell_planarise`),
`callback(volmesh, k, callback_args)` (`volmesh_planarise`), `callback(volmesh, formdiagram, k, callback_args)`
(`volmesh_reciprocate`) and `callback(cell, callback_args)` (`cell_arearise_face`): record frames there
to animate with `@viewer.on`. Examples: `examples/cell_manipulation/`, `examples/multi_cell_polyhedron/`.
