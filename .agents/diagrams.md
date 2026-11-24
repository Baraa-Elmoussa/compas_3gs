# Diagrams and datastructures

`from compas_3gs.diagrams import Cell, EGI, FormNetwork, ForceVolMesh, CellNetwork`
`from compas_3gs.datastructures import Mesh3gs, Network3gs, VolMesh3gs`

Every class is a compas 2 datastructure with all of its methods; the table lists only what
`compas_3gs` adds. Constructors are the compas ones (`from_vertices_and_faces`,
`from_vertices_and_cells`, `from_lines`, `from_meshgrid`, ...).

| Class | Base | Is | Adds |
|---|---|---|---|
| `Mesh3gs` | `compas.datastructures.Mesh` | base mesh | `face_area(f)` (vector area), `face_normal(f, unitized=True)` (a `Vector`, with a fallback for degenerate faces), `vertex_update_xyz(v, xyz, constrained=True)` (respects vertex attributes `x_fix y_fix z_fix`) |
| `Cell` | `Mesh3gs` | one polyhedral cell: the force polyhedron of one node | `cell_center()` (volume centroid); face attributes `target_area`, `target_normal` |
| `EGI` | `Mesh3gs` | extended Gaussian image of a cell: one vertex per face normal on the unit sphere | `EGI.from_volmesh_cell(ckey, volmesh)`; attribute `origin`; vertex attributes `normal`, `type` |
| `Network3gs` | `compas.datastructures.Graph` | base network | `vertex_update_xyz(node, xyz, constrained=True)`, `edge_avg_length()` |
| `FormNetwork` | `Network3gs` | form diagram | supports, forces, results: see [`model.md`](model.md); `dual` (its `ForceVolMesh`), `dual_halfface(edge)`, `update_angle_deviations()`; edge attributes `l_min l_max target_vector target_length` (used by reciprocation) |
| `VolMesh3gs` | `compas.datastructures.VolMesh` | base volmesh | `cell_pair_halffaces(c1, c2)` (the two halffaces where two cells touch, `KeyError` if they do not), `vertex_normal(v)`, `vertex_update_xyz(v, xyz, constrained=True)` |
| `ForceVolMesh` | `VolMesh3gs` | polyhedral force diagram (multi-cell) | `primal` (its `FormNetwork`); attributes `scale`, `convention` (-1: cells outward = members in compression); face attributes `target_area`, `target_normal`, `_a` |
| `CellNetwork` | `Network3gs` | disjointed force diagram: a cell per node (node attribute `cell`) | nothing builds one yet; `GSModel.cells` holds the cells of a form-first model |

## Correspondence form ↔ force

- A cell of `ForceVolMesh` is a node of `FormNetwork` (same key, after `volmesh_dual_network`).
- An interior face (a pair of halffaces) is a member: the member is perpendicular to the face and
  its force is the face area. Boundary halffaces are the external forces.
- `utilities.pair_hf_to_uv(volmesh, network)` / `pair_uv_to_hf(network, volmesh)` map halffaces and edges.
- `operations.point_reflection(volmesh)` turns compression into tension: the cells turn inside out
  while every face normal keeps its direction.

## Sample data (`compas_3gs.DATA`, the repository `data/` folder)

| File | Format | Build with |
|---|---|---|
| `cell_gss.json` | `{"vertices", "faces"}` | `Cell.from_vertices_and_faces(d["vertices"], d["faces"])` |
| `cells.json` | `{"cells": [{"vertices", "faces"}, ...], "point"}` (5 cells) | the same, per cell |
| `equilibrated_force_vectors.json` | `[{"origin", "vectors"}, ...]` (8 groups in equilibrium) | `egi_from_vectors(dict(enumerate(vectors)), origin)` |
| `volmesh_*.json` (`cubes_4x4`, `donut`, `honeycomb`, `mat`, `shells`, `warped`) | `{"vertices", "cells"}` | `ForceVolMesh.from_vertices_and_cells(d["vertices"], d["cells"])` |

Volmeshes are stored as plain vertices and cells, not with `to_json`: see [`compas2.md`](compas2.md).
