# Operations

`from compas_3gs.operations import ...` — topological and geometric changes of cells and volmeshes.
They change their input in place; most also return it.

## Cells

| Function | Does |
|---|---|
| `cell_split_indet_face_vertices(cell, fkey)` | split the vertices of a face that have more than one trailing edge, so the face can be pulled without tilting its neighbours |
| `cell_relocate_face(cell, fkey, xyz, normal)` | move a face to a target plane: each of its vertices slides along its trailing edge, so the neighbouring faces keep their planes (indeterminate vertices are split first) |
| `cell_collapse_short_edge(cell, u, v, min_length=0.1)` | move both ends of a short edge to its midpoint (the topology stays) |
| `cell_face_subdivide_barycentric(cell, fkey)` | fan of triangles around the face centroid |
| `cell_merge_coplanar_adjacent_faces(cell, tol=None)` | merge neighbouring coplanar faces |

## Volmeshes

| Function | Does |
|---|---|
| `volmesh_vertex_lift(volmesh, vkey, target_xyz, hfkeys)` | duplicate a boundary vertex and add cells between the given halffaces and the copy |
| `volmesh_vertex_merge(volmesh, vkeys, target_xyz)` | merge vertices into one at a point |
| `volmesh_halfface_pinch(volmesh, hfkey, xyz)` | add a pyramid cell on a boundary halfface |
| `volmesh_merge_adjacent_halffaces(volmesh, hfkeys)` | add a cell on a chain of adjacent boundary halffaces (their convex hull) |
| `volmesh_cell_subdivide_barycentric(volmesh, ckey)` | split a cell into pyramids around its centre; returns the new cell keys |
| `point_reflection(datastructure)` | reflect a graph, mesh or volmesh through the centroid of its vertices; on a force diagram this switches compression and tension (cells inside out, normals kept) |
| `volmesh_halfface_extrude`, `volmesh_halfface_subdivide`, `volmesh_cell_split_vertex` | raise `NotImplementedError` |

Each volmesh operation rebuilds the cells it changes through the public `VolMesh` API and copies
their cell, face and edge attributes back. Examples: `examples/cell_manipulation/00_10_cell_pull_face.py`,
`examples/multi_cell_polyhedron/01_60_volmesh_modification.py`, `examples/model/03_force_first.py`.
