# compas 2.15.1: what to use, and what is broken

`compas_3gs` uses compas 2 methods wherever they exist. These are the facts that are not visible
from the signatures; the upstream bugs are also listed in `TODO.md`.

## Broken in compas 2.15.1 — do not use

| Call | Problem | Use instead |
|---|---|---|
| `VolMesh` `to_json` / `from_json` / `copy()` | `__data__` writes each face once per halfedge: a 2×2×2 grid comes back with 192 halffaces instead of 48; raises `TypeError` once a cell attribute is set | store plain `{"vertices", "cells"}`, rebuild with `from_vertices_and_cells` |
| `volmesh.cell_vertex_faces`, `cell_vertex_neighbors` | `TypeError: 'dict_keys' object is not subscriptable` | `mesh = volmesh.cell_to_mesh(c)`: mesh vertex `j` is `cell_vertices(c)[j]`, mesh face `i` is `cell_faces(c)[i]` (see `EGI.from_volmesh_cell`) |
| `volmesh.halfface_manifold_neighbors`, `halfface_manifold_neighborhood` | `KeyError` on every boundary halfface | no replacement |
| `volmesh.is_valid()` | `NotImplementedError` | every `cell_to_mesh(c).is_closed()`; signed `volume_polyhedron(volmesh.cell_to_vertices_and_faces(c)) > 0`; halffaces with the same vertices are each other's opposite |
| `mesh.merge_faces([a, b, c, ...])` | wrong with three or more faces | merge two at a time; the merged face gets a new key |

## Surprising — use with care

| Call | Behaviour |
|---|---|
| `Mesh.volume(unify_cycles=False)` | the absolute volume: an inside-out cell is positive. For orientation use `volume_polyhedron(mesh.to_vertices_and_faces())` |
| `volmesh.is_halfface_on_boundary(hf)` | true when the halfface has no opposite cell, so it cannot check pairing |
| `VolMesh.edge_cells((u, v))` | the cells around an edge in dict order, not cyclic order; walk around the edge through `edge_halffaces` (see `volmesh_dual_volmesh`) |
| `VolMesh.edge_halffaces((u, v))` | only the halffaces running along the edge in that direction; look up `(v, u)` too |
| `VolMesh.cell_neighbors(c)` | lists a neighbour once per shared halfedge: de-duplicate with `dict.fromkeys(...)` |
| `VolMesh.delete_cell(c)` | also deletes the edge data of the cell's edges; when rebuilding a cell, copy cell, face and edge attributes back |
| `VolMesh.faces()` / `halffaces()` | one halfface per face / all halffaces |
| `VolMesh.face_area`, `halfface_area`, `Mesh3gs.face_area` | the vector area of a (possibly warped) face |
| `Mesh.centroid()` | area-weighted; for the centre of the vertices use `centroid_points(...)` |
| `mesh.face_neighbors(f)` | can list a face twice; de-duplicate before deleting or merging |
| `mesh.dual(cls=...)` | keeps keys but its faces point inward |
| `graph.edge_direction((u, v))` / `edge_vector((u, v))` | unit vector / full vector from `u` to `v` |
| `add_vertex(key=...)`, `add_cell(faces, ckey=...)` | the keyword is `key`; an unknown keyword (`vkey=`) is stored as an attribute and compas picks its own key |
| keys | not indices: after deletions they have gaps; `cell_to_mesh`, `to_vertices_and_faces`, `from_vertices_and_faces` renumber |

## Conventions in this package

- Tolerances from `compas.tolerance.TOL`. Two computed points are the same when closer than
  `TOL.angular` (unit sphere: EGI normals) or `TOL.absolute` (model coordinates), never by rounded
  geometric keys.
- Colours: `Color.from_i(i).rgb255` for gradients; 0–255 tuples or `Color` objects, never float tuples.
- No matplotlib and no Rhino: drawing is compas_viewer only ([`viewer.md`](viewer.md)).
