# compas 2.15 and 2.16: what to use, and what is broken

`compas_3gs` uses compas 2 methods wherever they exist. These are the facts that are not visible
from the signatures; the upstream bugs are also listed in `TODO.md`. compas 2.16.0 changes none of them.

## Broken in compas `VolMesh` — use a `VolMesh3gs`

`VolMesh3gs`, and so `ForceVolMesh`, replaces these methods. On a plain compas `VolMesh` they still fail.

| Call | On a compas `VolMesh` | On a `VolMesh3gs` |
|---|---|---|
| `to_json` / `from_json` / `copy()` | `__data__` writes each face once per halfedge: a 2×2×2 grid comes back with 192 halffaces instead of 48; raises `TypeError` once a cell attribute is set | keeps the vertex, halfface and cell keys and every attribute; also reads the data of a compas `VolMesh` |
| `cell_vertex_faces(c, v)`, `cell_vertex_neighbors(c, v)` | `TypeError: 'dict_keys' object is not subscriptable` | in order around the vertex; face `i` has the halfedge from the vertex to neighbour `i` |
| `halfface_manifold_neighbors`, `halfface_manifold_neighborhood` | `KeyError` on every boundary halfface | the halffaces of the neighbouring cells that continue the halfface across its edges |
| `is_valid()` | `NotImplementedError` | checks the topology: halffaces in cells, closed cells, matching neighbours. It does not check orientation: an inside-out cell has a negative `volume_polyhedron(volmesh.cell_to_vertices_and_faces(c))` |
| `edge_halffaces((u, v))`, `edge_cells((u, v))` | in dict order, not in order around the edge | in order around the edge, from the boundary for an edge on the boundary; still one halfface per cell, the one with the halfedge `(u, v)` |
| `delete_cell(c)` | leaves the face and edge attributes of the cell | removes the attributes of the faces and edges that no other cell uses; when rebuilding a cell, copy cell, face and edge attributes back |

Code that takes any compas `VolMesh` (`volmesh_dual_volmesh`, `EGI.from_volmesh_cell`) walks around the
edge itself and goes through `cell_to_mesh` instead.

## Broken in compas `Mesh` — do not use

| Call | Problem | Use instead |
|---|---|---|
| `mesh.merge_faces([a, b, c, ...])` | wrong with three or more faces | merge two at a time; the merged face gets a new key |

## Surprising — use with care

| Call | Behaviour |
|---|---|
| `Mesh.volume(unify_cycles=False)` | the absolute volume: an inside-out cell is positive. For orientation use `volume_polyhedron(mesh.to_vertices_and_faces())` |
| `volmesh.is_halfface_on_boundary(hf)` | true when the halfface has no opposite cell, so it cannot check pairing |
| `VolMesh.edge_halffaces((u, v))` | only the halffaces running along the edge in that direction; look up `(v, u)` too |
| `VolMesh.cell_edges(c)` | the halfedges of the cell, so every edge twice, once each way |
| `VolMesh.cell_neighbors(c)` | lists a neighbour once per shared halfedge: de-duplicate with `dict.fromkeys(...)` |
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
